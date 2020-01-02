#!/lab/corradin_data/FOR_AN/anaconda3/bin/python
#/usr/bin/env python

"""
Orchestrator takes the super_pipeline configuration and manages all the jobs
for running the analysis against a single chromosome.

* Takes the same arguments as the super_pipeline.


Usage:
	./orchestrator.py file1 file2 pairing_file init_file p_file [odds_file]

	The meanings of file1 and file2 change depending on which pipeline is
	configured in the init_file.

Overview:

* The orchestrator reads in the pairing_file and creates an overall task for
	each rsid pair.

* A task contains one or more jobs. A job contains the status of a single call
	out to the super_pipeline to do some number of randomizations for a single
	rsid_pair.

* The orchestrator loops (checking every minute) to see the status of jobs. If
	a job fails a new job is created to continue the work.

* The orchestrator's state information (Tasks & Jobs) are written to written to
	disk. So the orchestrator itself can recover if it crashes.

"""

from contextlib import contextmanager
from datetime import datetime
from types import SimpleNamespace
from collections import namedtuple
from atomicwrites import atomic_write

import subprocess
from subprocess import Popen, PIPE

import argparse
import glob
import logging
import math
import os
import pandas as pd
import pickle
import re
import shutil
#import sys
import time
#import uuid
import numpy as np
from functools import wraps
from toolz import partial
#from gooey import Gooey
from num2words import num2words


MAX_JOB_RETRY_COUNT = 10
MAX_ITERATIONS_PER_JOB = 100_000
RsidPair = namedtuple("RsidPair",["GWAS_rsid", "outside_rsid"])
SinglePairing = namedtuple("SinglePairing",["pairing_file_line", "rsid_pair_tuple"])

logging.basicConfig(level=logging.INFO)

@contextmanager
def cd(newdir):
	"""
	mimicking the cd command in terminal
	"""
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	#except FileNotFoundError:
		#raise FileNotFoundError(f"Does not exist path {newdir} from directory {prevdir}")
	finally:
		os.chdir(prevdir)

@contextmanager
def open_with_exception(file, *args, **kwargs):
	try:
		with open(file, *args, **kwargs) as f:
			yield f
	except FileNotFoundError:
		raise FileNotFoundError(f"Cannot find file {file}")

def checkpoint(func):
	@wraps(func)
	def wrapped(*args, **kwargs):
		instance = args[0]
		klass = instance.__class__
		func_ret = func(*args, **kwargs)

		klass.write_to_file(instance)
		return func_ret
	return wrapped


class RsidTask(SimpleNamespace):
	"""
	Represents status of ALL the work that needs to be done for a single rsid_pair.
	* Has many jobs.

	statuses:
		new -> running -> complete
	"""
	# rsid_pair: tuple
	# pipeline_arguments: list
	# target_iterations: int
	# status: str
	# # completed_iterations: int
	# jobs: list  # of JobStatus objects

	def stop_monitoring(self, job):
		# Only the jobs list is actively monitored.
		# Moving the job into the ended_jobs list stops the monitoring.
		self.ended_jobs.append(job)
		self.jobs.remove(job)

	def rendered_command(self, partition, iterations):
		"""
		Returns fully specified pipeline command with all the arguments and paths.
		E.g. TODO: Add example
		:param partition: Which batch of the jobs. This translates into the filesystem path /outside/inside/<batch>/job/*.*
		:param iterations: How many iterations is this job doing?
		:return:
		"""
		#GWAS_rsid, outside_rsid = self.rsid_pair
		rsid_pair = self.rsid_pair
		rendered_command = self.pipeline_script_with_args.format(GWAS_rsid=rsid_pair.GWAS_rsid,
													 outside_rsid=rsid_pair.outside_rsid,
													 pairing_file_line=self.single_pairing.pairing_file_line,
													 target_num_iterations=self.target_iter_str,
													 job_num_iterations=iterations,
													 partition=partition)

		return rendered_command

	def rendered_job_name(self, partition, iterations):
		"""
		Returns the lsf job name
		e.g. TODO: Add example
		:param partition:
		:param iterations:
		:return:
		"""
		# should I do this or just split at "/" (get the final part of output_path)?
		#GWAS_rsid, outside_rsid = self.rsid_pair
		rsid_pair = self.rsid_pair

		return self.job_name_template.format(GWAS_rsid=rsid_pair.GWAS_rsid,
											 outside_rsid=rsid_pair.outside_rsid,
											 target_num_iterations=self.target_iter_str,
											 job_num_iterations=iterations,
											 partition=partition)

	# set the folder_path to descriptor object so that the folder path can change
	# every time the task's attribute changes without explicitly changing it.
	# Useful when updating task's target iteration
	# can make target_iterations a descriptor instead and change folder path in __set__() method

	# dynamically fetches task's attribute to prevent stale data when one attribute changes:

	@property
	def rsid_pair(self):
		return self.single_pairing.rsid_pair_tuple

	@property
	def name(self):
		rsid_pair = self.rsid_pair
		return self.task_folder_template.format(GWAS_rsid=rsid_pair.GWAS_rsid,
												outside_rsid=rsid_pair.outside_rsid,
												target_num_iterations=self.target_iter_str)
	@property
	def folder_path(self):
		formatted_task_path = "{}/{}".format(self.rsid_pair_path, self.name)
		return formatted_task_path

	@property
	def rsid_pair_path(self):
		rsid_pair = self.rsid_pair
		return self.outside_folder_struct_template.format(GWAS_rsid=rsid_pair.GWAS_rsid,
											  outside_rsid = rsid_pair.outside_rsid)
	@property
	def target_iter_str(self):
		return num2words(self.target_iterations).replace(" ", "_")



"""
The status of a single process (e.g. one LSF Job)

statuses:
	new -> running -> complete
				  \-> failed -> failed_requeued
"""
class Job(SimpleNamespace):
	def stub():
		"""
		This is just so we can subclass SimpleNamespace.
		Do don't have any custom methods needed yet.
		"""
		return False
#    rsids: tuple
#    command: str
#    jobid: int
#    status: str
#    # audit columns? created_at, started_at, completed_at?

class LsfJobRunner:
	def __init__(self, **kwargs):
		valid_attributes = ["queue_name", "job_group", "outside_folder_struct_template"]
		for attr_name, attr_val in kwargs.items():
			if attr_name in valid_attributes:
				setattr(self, attr_name, attr_val)
			else:
				raise ValueError("Unrecognized attribute: ", attr_name)
		# self.queue_name = kwargs["queue_name"]
		# self.job

	def run(self, job):
		"""
		Queue this job in LSF.
		:param job: The Job object containing the command we want to run.
		:return: job
		"""
		if job.command is None:
			raise ValueError("Job command cannot be None. Please check rendered_command method or your input")

		filled_outside_folder_struct = self.outside_folder_struct_template.format(GWAS_rsid=job.rsid_pair[0], outside_rsid=job.rsid_pair[1])
		# save template for debugging using email
		#-o {filled_outside_folder_struct}/myStdOut.out -e {filled_outside_folder_struct}/myStdErr.err
		cmd = "bsub -o /dev/null -e {filled_outside_folder_struct}/myStdErr.err -q {queue_name} -g {job_group} -J {job_name} {job_command}".format(filled_outside_folder_struct=filled_outside_folder_struct,queue_name = self.queue_name, job_group=self.job_group, job_name=job.name, rsid_pair=job.rsid_pair, job_command=job.command)

		logging.info("RUN: %s" % cmd)
		bsub_output = subprocess.check_output(cmd, shell=True)
		print(bsub_output)
		job.jobid = self.__parse_jobid(str(bsub_output))
		job.status = "running"
		return job

	def check_job_status(self, job_folder, jobid):
		"""
		Determine the job status from LSF.
		:param jobid:
		:return: "running" | "complete" | "failed" | "unknown"
		"""
		# try:
		bjobs_cmd = ["bjobs", jobid]
		logging.info("\tbjobs command: '%s'" % bjobs_cmd)
		p = Popen(bjobs_cmd, stdout=PIPE, stderr=PIPE)
		bjobs_stdout, bjobs_stderr = p.communicate()
		#bjobs_stdout = (stdout)
		#bjobs_stderr = str(stderr)

		logging.debug("\tbjobs stdout: %s" % (bjobs_stdout))
		logging.debug("\tbjobs stderr: %s" % (bjobs_stderr))
		if re.compile("Job <\\d+> is not found").search(str(bjobs_stderr)) is None:
			# Found in bjobs output
			#lines = str(bjobs_stdout).split("\\n")
			lines = bjobs_stdout.splitlines()
			logging.info("lines: %s" %  lines)
			logging.info("len lines: %d" %  len(lines))
			if len(lines)>=2:
				# JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
				# 5612457 andy    RUN   normal     tak4        it-c03b05   sleep 1000 Jan  8 15:46
				statline = str(lines[1])
				statline_arr = re.compile("\s+").split(statline)
				logging.debug(statline_arr)
				stat = statline_arr[2]
				logging.debug("line1: %s" % lines[1])
				logging.debug("stat: %s" % stat)
				if stat == "PEND":
					return "running"
				elif stat == "RUN":
					return "running"
				elif stat == "DONE":
					return "complete"
				elif stat == "EXIT":
					return "failed"
		else:
			# # Not found in bjobs output
			# # Either too old or invalid id
			# # check bhist
			# logging.debug("\tNot found with bjobs")

			# bhist_cmd = ["bhist", "-b", jobid]
			# logging.debug("\tbhist command: '%s'" % bhist_cmd)
			# p = Popen(bhist_cmd, stdout=PIPE, stderr=PIPE)
			# stdout, stderr = p.communicate()
			# bhist_stdout = str(stdout)
			# bhist_stderr = str(stderr)
			# logging.debug("\tbhist stdout: %s" % (bhist_stdout))
			# logging.debug("\tbhist stderr: %s" % (bhist_stderr))

			# if re.compile("Done successfully\.").search(bhist_stdout):
			# 	return "complete"

			logging.debug("\tCannot find job in bjobs checking for \'finished\' file in job folder")
			logging.debug("\tThe current folder is {}".format(os.getcwd()))
			logging.debug("\tThe job folder is {}".format(job_folder))

			try:
				with cd(job_folder):
					if os.path.isfile("finished"):
						logging.debug("\tfound file \'finished\'")
						return "complete"
					else:
						logging.debug("\tcannot find file \'finished\' in directory {}. Marking job as failed".format(
							os.getcwd()))
						return "failed"

			except FileNotFoundError:
				logging.info("\tin except clause b/c FileNotFoundError exception in check_job_status_function, returning failed".format(
					os.getcwd()))
				return "failed"
		# except FileNotFoundError:
		# 	# Probably no bjobs.
		# 	logging.warn("\tFileNotFound when checking status. Probably no bjobs command. Returning unknown")
		# 	return "unknown"
		# logging.warn("\tUnhandled bjobs output. Returning unknown.")
		return "unknown"  # parse bstatus

	def __parse_jobid(self, bsub_output):
		jobid = re.search("<(\d+)>", bsub_output).group(1)
		return jobid

	def cleanup_job(self, job):
		cmd = f"bkill {job.jobid}"
		print(cmd)
		bsub_output = subprocess.check_output(cmd, shell=True)
		print(bsub_output)

# class LocalJobRunner:
#     def __init__(self, **kwargs):
#         self.threads_count = kwargs["threads_count"]
#         self.job_queue = []
#         self.start_workers()
#
#     def run(self, job):
#         self.job_queue.append(job)
#         self.jobid = uuid.uuid4()
#         self.status = "running"
#         return job
#
#     def check_job_status(self, jobid):
#         return "running"
#
#     def start_workers(self):
#         # spin up n threads (self.threads_count)
#         # multiprocess Process, Pool etc.
#         return True

# TODO: Talks more about task than self, move into Task

def create_new_job(task, partition, job_iterations, status):
	"""
	Creates a Job record from a task.
	partition is for filesystem organization (i.e. a batch of jobs)
	iteration is the number of iterations this job will perform.
	"""
	job = Job(rsid_pair=task.rsid_pair,
			   command=task.rendered_command(partition, job_iterations),
			   name=task.rendered_job_name(partition, job_iterations),
			   fail_count=0,
			   status=status)

	job.job_folder_path = "{}/{}_all_files".format(task.folder_path, job.name)
	return job

# TODO: Talks more about task than self, move into Task
def create_new_jobs(task,status):
	"""
	Creates one or more Job to accomplish at least task.target_iterations
	"""
	iterations_per_job = task.target_iterations
	if iterations_per_job > MAX_ITERATIONS_PER_JOB:
		iterations_per_job = MAX_ITERATIONS_PER_JOB
		jobs_count = math.ceil(task.target_iterations / iterations_per_job) # If there is a remainder, run an extra job to cover.
		jobs = [create_new_job(task, i, iterations_per_job, status) for i in range(jobs_count)]

	#if the target iteration is smaller than MAX_ITERATIONS_PER_JOB then don't split and indicate in job name that there is only one partition
	else:
		jobs = [create_new_job(task,"all", iterations_per_job, status)]


	task.jobs.extend(jobs)
	return jobs

def read_pairs(pair_file, cli_args):
	"""
	:param pair_file: - Same file taken by the super_pipeline script. One-line per pair, space separated.
	:return: List of tuples of the pairs.
	"""
	single_pairings_list = []
	rsid_pairs_list = []
	with open_with_exception(pair_file, 'r') as pairing_file:
		for line in pairing_file:
			clean_line = line.strip()
			line_list = clean_line.split(" ")
			#CANNOT HANDLE TAB

			if cli_args.pipeline_type == "TRANS_CC":
				rsid_pair = RsidPair(line_list[0], line_list[3])
			else:
				rsid_pair = RsidPair(*line_list)

			single_pairing = SinglePairing(clean_line, rsid_pair)
			rsid_pairs_list.append(rsid_pair)
			single_pairings_list.append(single_pairing)

		if len(rsid_pairs_list) != len(set(rsid_pairs_list)):
			raise ValueError("The pairing file has duplicated rsid pairs, please remove the duplicate(s) and try again")
	return single_pairings_list

class JobOrchestrator:
	tasks_pickle_file = "orchestrator_tasks.pickle"
	orchestrator_pickle_file = "orchestrator.pickle"
	z_dict = {99: 2.576, 98: 2.326, 95: 1.96, 90: 1.645, 80: 1.28}

	@classmethod
	def load_from_file(cls):
		print(f"Loading orchestrator from file {cls.orchestrator_pickle_file}")
		with open_with_exception(cls.orchestrator_pickle_file, 'rb') as f:
			orchestrator = pickle.load(f)
			assert orchestrator.__class__ == cls
			return orchestrator

	@classmethod
	def write_to_file(cls, instance):
		print(f"Checkpointing orchestrator to file {cls.orchestrator_pickle_file} in directory {os.getcwd()}")
		with atomic_write(cls.orchestrator_pickle_file, mode='wb+', overwrite=True) as f:
			pickle.dump(instance, f)


	@checkpoint
	def __init__(self, z_threshold, single_pairings, pipeline_args, outside_folder_struct_template, task_folder_template, job_folder_template, args):
		#addded args
		self.single_pairings = single_pairings
		self.pipeline_args = pipeline_args
		self.tasks = {}
		self.completed_tasks = {}
		self.z = self.__class__.z_dict[int(z_threshold)]
		# placing the runner creation here so can pass in user input easily: queue name and job group
		self.runner = LsfJobRunner(queue_name=args.queue, job_group=args.group, outside_folder_struct_template=outside_folder_struct_template)
		#end

		self.overall_start_time = time.time()

		# Total number of iterations to do at each filter_step.

		self.iteration_steps = list(np.logspace(3,args.max_iterations,num=args.max_iterations -3 + 1).astype(int))
		print(self.iteration_steps)

		#10_000_000# TODO: param
		self.highest_iteration = self.iteration_steps[0]
		# Also see MAX_ITERATIONS_PER_JOB

		pipeline_script = "/lab/corradin_data/FOR_AN/OUTSIDE_VARIANT_PIPELINE/github_repos/outside-variants/python_scripts/OVP_script.py -op {GWAS_rsid}_{outside_rsid}"
		pipeline_script_with_args = "{} {}".format(pipeline_script, " ".join(pipeline_args))

		with open("pipeline_info", "w+") as f:
			f.write(pipeline_script_with_args)
			f.write('\n')
			f.write('\n')
			f.write("Iteration steps: {}".format(str(self.iteration_steps)))
			f.write('\n')
			f.write('\n')
			f.write("MAX_JOB_RETRY_COUNT = {}, MAX_ITERATIONS_PER_JOB = {}".format(MAX_JOB_RETRY_COUNT, MAX_ITERATIONS_PER_JOB))
			print("Pipeline script with args: ", pipeline_script_with_args)

		# initialize Tasks, make sure we have one per rsid_pair.
		# TODO: make sure the target iteration is correct when orchestrator is loaded and not init from scratch
		for single_pairing in single_pairings:
			rsid_pair = single_pairing.rsid_pair_tuple
			#create folder structure so error files for each pair can be recorded
			os.makedirs(outside_folder_struct_template.format(GWAS_rsid=rsid_pair.GWAS_rsid, outside_rsid=rsid_pair.outside_rsid), exist_ok=True)
			task = self.tasks.get(rsid_pair, None)
			if task is None:
				task = RsidTask(single_pairing = single_pairing,
								pipeline_script_with_args=pipeline_script_with_args,
								outside_folder_struct_template= outside_folder_struct_template,
								task_folder_template = task_folder_template,
								job_name_template = job_folder_template,
								target_iterations = self.iteration_steps[0],
								status="new",
								jobs=[],
								ended_jobs=[])
			self.tasks[rsid_pair] = task
			self.completed_tasks[rsid_pair] = []


		#self.__class__.write_to_file(self)

	def run(self):
		"""
		Loop until all the rsid_pairs are completely processed.
		Starting & restarting jobs as necessary.
		"""

		i = 0
		with open("job_log", "a+") as job_log, open("task_log", "a+") as task_log:
			labels = ["GWAS_rsid", "outside_rsid", "task_name", "job_name", "status"]
			job_log.write("\t".join(labels))
			job_log.write("\n")

			task_log.write("\t".join(labels))
			task_log.write("\n")

			while self.incomplete(self.tasks):
				done_tasks = []
				print(f"Checked {i} times")
				i +=1

				for rsid_pair in self.tasks:
					task = self.tasks.get(rsid_pair, None)
					logging.info("rsid_pair %s,%s" % rsid_pair)

					# First run initialization of jobs.
					if len(task.jobs) == 0:
						logging.info("\tstarting first job")
						new_jobs = create_new_jobs(task, "new")
						for job in new_jobs:
							self.runner.run(job)
						task.status = "running"

					# Re-check all the jobs for the task.

					task.all_done = self.check_task_jobs(job_log=job_log, task= task)

					# Split child jobs
					if task.all_done:

						line = [f"{task.rsid_pair.GWAS_rsid}",f"{task.rsid_pair.outside_rsid}",f"{task.name}", "NA"]
						task.need_split_cleaned_up = self.needs_split(task)
						if task.need_split_cleaned_up:
							current_index = self.iteration_steps.index(task.target_iterations)
							if current_index+1 > len(self.iteration_steps) - 1:
								logging.info("MAX ITERATION REACHED, STILL NEED MORE PERM FOR RSID PAIR {} AT {} ITERATIONS".format(task.rsid_pair,task.target_iter_str))
								# remove task and move on to next task
								line.append("reached_max_iter_more_perm")
								task_log.write("\t".join(line))
								task_log.write("\n")
								done_tasks.append(task)
								continue
							else:
								# try to move to the next iteration step
								task.target_iterations = self.iteration_steps[current_index + 1]
								logging.info(
									f"MOVING TO NEXT STEP OF {task.target_iter_str} ITERATIONS, STILL NEED MORE PERM FOR RSID PAIR {task.rsid_pair} AT {num2words(self.iteration_steps[current_index])} ITERATIONS")

								#update highest iteration:
								if task.target_iterations > self.highest_iteration:
									self.highest_iteration = task.target_iterations

								#create new jobs and run them
								next_iter_step_jobs = create_new_jobs(task, "new")
								for job in next_iter_step_jobs:
									self.runner.run(job)
						else:
							logging.info("DONE WITH RSID PAIR {} AT {} ITERATIONS".format(task.rsid_pair, task.target_iter_str))
							task.status = "complete"
							line.append(f"complete_{task.target_iter_str}")
							task_log.write("\t".join(line))
							task_log.write("\n")
							#self.stop_monitoring(task)
							done_tasks.append(task)

					print("-")
				print("---")
				# print(self.tasks)
				print("===")
				logging.info(f"Currently in this directory: {os.getcwd()}")

				#removing all the done tasks at once:
				for finished_task in done_tasks:
					checkpoint(self.stop_monitoring(finished_task))
				#self.save_tasks()
				time.sleep(60)

		self.final_combine()
		print("all done ---------------")
		self.overall_end_time = time.time()
		print(f"Finished {len(self.single_pairings)} SNP pairs from {self.iteration_steps[0]} to {self.highest_iteration} in {self.overall_end_time - self.overall_start_time}")


	@checkpoint
	def check_task_jobs(self, *, job_log, task):
		all_done = True
		for job in task.jobs:
			logging.info(job.name)

			# Recheck job status

			logging.debug("\tjob status was %s" % (job.status))
			job.status = self.runner.check_job_status(job.job_folder_path, job.jobid)
			logging.debug("\tnew status is %s" % (job.status))


			if job.status == "running":
				logging.info("\tJob still running")
				all_done = False

			elif job.status == "complete":
				logging.info("\tJob complete!")
				task.stop_monitoring(job)

			else:
				all_done = False
				line = [f"{job.rsid_pair[0]}", f"{job.rsid_pair[1]}", f"{task.name}", f"{job.name}"]
				if job.status == "failed":
					logging.info("\tJob failed! Restarting")
					job.fail_count += 1

					if job.fail_count > MAX_JOB_RETRY_COUNT:
						logging.warn("\tExceeded Maximum Retry count for. Giving up!")
						# TODO: How do we know this happened without looking at the logs?
						line.append("failed_max_retry_reached")
						task.stop_monitoring(job)
					else:
						self.runner.run(job)
					line.append("failed_restarted")

				else:
					logging.warn("unhandled status " + job.status)
					line.append("unhandled_status")
					raise RuntimeError("Unknown status encountered")

				job_log.write("\t".join(line))
				job_log.write("\n")
		return all_done

	def stop_monitoring(self, task):
		# Only the tasks list is monitored.
		self.completed_tasks[task.rsid_pair] = task
		del self.tasks[task.rsid_pair]

	def incomplete(self, tasks):
		num_incomplete = 0
		for task in tasks:
			if task != "complete":
				num_incomplete +=1

		if num_incomplete:
			print("-----")
			print("----")
			print("---")

			print(f"Still have {num_incomplete} tasks incomplete")
			return True
		else:
			return False

	# def save_tasks(self):
	# 	with open(self.tasks_pickle_file, 'wb') as f:
	# 		pickle.dump(self.tasks, f)

	# def load_tasks(self):
	# 	with open_with_exception(self.tasks_pickle_file, 'rb') as f:
	# 		return pickle.load(f)


	def needs_split(self,task):  # STUB
		"""
		Inspects the checkpoint files to determine if the task is inconclusive and needs to be split and have more iterations run.
		:param task:
		:return: True | False
		"""
		# read checkpoint files
		# math to determine confidence
		def compare_mtc_threshold(mtc_threshold, combined_iter_used_for_p_value, combined_iterations):
			combined_p_value = combined_iter_used_for_p_value/combined_iterations
			p_value_no_zero = combined_p_value if combined_p_value != 0 else 1/combined_iterations
			factor = self.z * np.sqrt((p_value_no_zero * (1 - p_value_no_zero)) / combined_iterations)
			ci_neg, ci_pos = p_value_no_zero - factor, p_value_no_zero + factor

			if ((ci_neg == np.float64(0)) & (ci_pos == np.float64(0))) or ci_neg < mtc_threshold < ci_pos:
				status = "more_perm"

			elif mtc_threshold < ci_neg:
				status = "non_sig"

			elif mtc_threshold > ci_pos:
				status = "sig"

			else:
				status = "unknown"

			return  p_value_no_zero, ci_neg, ci_pos, status

		def remove_folder(path):
			# Try to remove tree; if failed show an error using try...except on screen
			try:
				shutil.rmtree(path)
			except OSError as e:
				print("Error: %s - %s." % (e.filename, e.strerror))

		#get all the finished jobs/partitions of a task
		with cd(task.folder_path):
			print("----")
			print(f"In task folder: {task.folder_path}")
			all_jobs_files = glob.glob(os.path.join('', "*_all_files/*_ALL_INFO*"))

			task_df = pd.concat(pd.read_csv(file, sep="\t") for file in all_jobs_files)

		with cd(task.rsid_pair_path):
			task_df.to_csv("{}_all_partitions".format(task.name), sep="\t",index = False)
			#ANDY_REVIEW: is writing a file and immediately reading it back in a good idea?
			all_tasks_so_far = glob.glob(os.path.join('', "*_all_partitions"))
			df = pd.concat(pd.read_csv(file, sep="\t") for file in all_tasks_so_far)

		# #find out which side we use to calculate p value, perm higher or perm lower
		# df["iter_used_for_p_value"] = np.where(df["odds_ratio_GWAS_and_outside"] > df["odds_ratio_GWAS"],
		#                                        df["num_perm_higher"], df["num_perm_lower"])

		# agg_df = df.groupby("GWAS_rsID,outside_rsID,GWAS_geno,outside_geno,mtc_threshold".split(","), as_index=True)

		# #aggregate stat from all jobs
		# agg_df = agg_df.agg({"iterations": "sum",
		#                      "iter_used_for_p_value": "sum",
		#                      "p_value_no_zero": "mean",
		#                      "CI_lower": "mean",
		#                      "CI_higher": "mean"}).add_prefix("combined_").reset_index()

		# filter_func = np.vectorize(compare_mtc_threshold)
		# agg_df["status"] = filter_func(agg_df["mtc_threshold"], agg_df["combined_CI_lower"],
		#                                agg_df["combined_CI_higher"])


		agg_df = df.groupby("GWAS_rsID,outside_rsID,GWAS_geno,outside_geno,mtc_threshold".split(","), as_index=True)
		agg_df = agg_df.agg({"iterations": "sum",
					"iter_used_for_p_value": "sum",
					 "num_perm_lower":"sum",
					 "num_perm_higher":"sum",
					}).add_prefix("combined_").reset_index()

		filter_func = np.vectorize(compare_mtc_threshold)
		#agg_df["status"],  = filter_func(agg_df["mtc_threshold"], agg_df["combined_iter_used_for_p_value"], agg_df["combined_iterations"])
		agg_df["combined_p_value_no_zero"], agg_df["combined_ci_neg"], agg_df["combined_ci_pos"], agg_df["status"] = filter_func(agg_df["mtc_threshold"], agg_df["combined_iter_used_for_p_value"], agg_df["combined_num_perm_lower"] + agg_df["combined_num_perm_higher"])

		info_cols = ['GWAS_rsID', 'outside_rsID', 'GWAS_geno', 'outside_geno', 'case_count',
			   'case_total', 'case_odds', 'control_count', 'control_total',
			   'control_odds', 'odds_ratio_GWAS', 'odds_ratio_GWAS_and_outside','outside_cause_higher_or_lower_risk']
		info_df = df[info_cols]
		agg_df = pd.merge(agg_df, info_df, on="GWAS_rsID,outside_rsID,GWAS_geno,outside_geno".split(","), how="left")

		agg_df = agg_df.sort_values(['GWAS_rsID', 'outside_rsID', 'GWAS_geno', 'outside_geno','combined_iterations'], ascending=False).drop_duplicates().sort_index()

		with cd(task.rsid_pair_path):
		# 	# clean up task folder
		# 	# TODO: Extract into cleanup method.
			if len(task.name) > 0:
				for old_task_folder in glob.glob( "*{}*/".format(task.name)):
					remove_folder(old_task_folder)

			#output combined data so far
			agg_df.to_csv("{}_combined".format(task.name), sep="\t", index = False)
			agg_df.to_csv("{}_{}_final".format(task.rsid_pair[0], task.rsid_pair[1]), sep="\t", index = False)
		return "more_perm" in agg_df["status"].unique()

	def final_combine(self):
		all_rsid_files = glob.glob(os.path.join('', "rs*/rs*/rs*_rs*_final"))
		all_rsid_df = pd.concat(pd.read_csv(file, sep="\t") for file in all_rsid_files)
		all_rsid_df = all_rsid_df.sort_values(['combined_iterations'], ascending=False).drop_duplicates().sort_values(['GWAS_rsID', 'outside_rsID', 'GWAS_geno', 'outside_geno'])
		all_rsid_df.to_csv("all_results", sep="\t", index = False)

	def cleanup_jobs(self):
		for rsid_pair in self.tasks:
			task = self.tasks.get(rsid_pair, None)
			for job in task.jobs:
				self.runner.cleanup_job(job)

	def restart_jobs(self):
		for rsid_pair in self.tasks:
			task = self.tasks.get(rsid_pair, None)
			for job in task.jobs:
				self.runner.run(job)


def run_setup(args):
	pairing_file = args.input_folder_path + "/" + args.SNP_pairs
	single_pairings = read_pairs(pairing_file, args)

	#print(args.to_list())
	if args.unique_identifier == args.input_folder_path.split("/")[-1]:
		raise ValueError("unique_identifier cannot be the same as input folder")

	outside_folder_struct_template = "{GWAS_rsid}/{outside_rsid}"
	task_folder_path_template = "{GWAS_rsid}_{outside_rsid}_" + args.unique_identifier + "_{target_num_iterations}_TargetIter"
	job_folder_path_template = task_folder_path_template + "_partition_{partition}"

	output_file_path_template = "{}/{}/{}".format(outside_folder_struct_template,task_folder_path_template, job_folder_path_template)

	pipeline_args = [args.input_folder_path, args.case_gen, args.control_gen, "\'{pairing_file_line}\'", args.init_file, output_file_path_template, args.unique_identifier, "--iter {job_num_iterations}"]

	if args.override:
		pipeline_args.append("--override")

	# threshold to calculate confidence intervals
	z_threshold = 99

	# orchestrator.submit_jobs(args.group)

	time_now = datetime.now().strftime("%Y-%m-%d")
	run_folder_name = f"{args.unique_identifier}_MAX_ITER_PER_JOB_{MAX_ITERATIONS_PER_JOB}_{time_now}"

	try:
		os.makedirs(run_folder_name, exist_ok = False)
	except FileExistsError:
		raise FileExistsError(f"Found folder of the same name: {run_folder_name}, either choose another name or delete this folder. If you want to rerun, change to rerun mode")

	with cd(run_folder_name):
		#saving the run command for future reference
		with open(f"run_{args.unique_identifier}.sh", "w+") as f:
		   f.write("#!/bin/bash")
		   f.write("\n")
		   f.write(" ".join(sys.argv))
		   f.write("\n")

		orchestrator = JobOrchestrator(z_threshold, single_pairings, pipeline_args, outside_folder_struct_template,
									   task_folder_path_template, job_folder_path_template, args)

	return run_folder_name, orchestrator


def rerun_setup(rerun_folder):
	with cd(rerun_folder):
		orchestrator = JobOrchestrator.load_from_file()
		orchestrator.iteration_steps = list(orchestrator.iteration_steps)
		return rerun_folder, orchestrator


def main():

	parser = argparse.ArgumentParser(description='orchestrator for outside var pipeline')
	subparsers = parser.add_subparsers(help="run modes", dest="run_mode")

	#rerun args
	parser_rerun = subparsers.add_parser('rerun', help='activate rerun mode')
	parser_rerun.add_argument('rerun_folder', help='folder to rerun')

	#run args
	parser_run = subparsers.add_parser('run_fresh', help='activate run mode')

	#get results
	parser_get_results = subparsers.add_parser('get_results', help='activate get results mode')
	parser_get_results.add_argument('rerun_folder', help='folder to rerun')

	parser_run.add_argument('input_folder_path', help='absolute path of input folder')
	parser_run.add_argument('case_gen', help='name chromosome-specific case .gen file, with path specified using argument PATH in init file. Example: ALL_MS_impute2_chr20.gen ')
	parser_run.add_argument('control_gen', help='chromosome-specific control .gen file, with path specified using argument PATH in init_file. Example: ALL_controls_58C_NBS_WTC2_impute2_chr20.gen')
	parser_run.add_argument('SNP_pairs', help='file with each line having GWAS_rsid and outside_rsid, separated by a delimiter. Example: SNPpairs_SAMPLE ')
	parser_run.add_argument('init_file', help='file with init arguments. Formatted in the form: keyword,value separated by tab')
	parser_run.add_argument('unique_identifier', help='unique identifier for this run. Will be used for creating names of folders and files')
	parser_run.add_argument('--override', action='store_true', help='if this argument is present, override any current folder with the same unique_identifier')
	parser_run.add_argument('--group', '-g', required=True, help= "job group to run this batch in. Use bgadd /[job_group_name] to add a job group. E.g: bgadd /MSlogreg" )
	parser_run.add_argument('--pipeline_type', '-pt',  nargs='?', const='other', default='other', choices=["TRANS_CC", "other"], help="specify type of pipeline, needed if run TRANS_CC pipeline since the pairing file format is different")
	parser_run.add_argument('--queue', nargs='?', const='all_corradin', default="all_corradin", choices=["all_corradin", "corradin", "normal"], help="which LSF queue to spawn jobs in")
	parser_run.add_argument('--max_iterations', nargs='?',type=int, default=7)


	args = parser.parse_args()
	# rerun_args = parser_rerun.parse_args()
	# run_args = parser_run.parse_args()
	# print(f"args: {args}, rerun_args: {rerun_args}, run_args: {run_args}")

	#testing args
	# args = parser.parse_args("/lab/corradin_data/FOR_AN/OUTSIDE_VARIANT_PIPELINE/for_pipeline_an/input_folder/MS ALL_MS_impute2_chr16.gen ALL_controls_58C_NBS_WTC2_impute2_chr16.gen Filtered_SNP_pairs_ALL_MS_impute2_chr16_ALL_controls_58C_NBS_WTC2_impute2_chr16_file Anna_init_get_mtc_logreg.txt logreg_100mil --override --group /test_orchestrator".split( ))

	if args.run_mode == "run_fresh":
		working_dir, orchestrator = run_setup(args)
	else:
		working_dir, orchestrator = rerun_setup(args.rerun_folder)

	if args.run_mode == "get_results":
		with cd(working_dir):
			orchestrator.final_combine()
	else:
		try:
			with cd(working_dir):
				orchestrator.run()
		except BaseException as e:
			print("\t ENCOUNTERED ERROR, cleaning up")
			#orchestrator.cleanup_jobs()
			print("exiting")
			raise e

if __name__ == '__main__':
	import sys
	main()




