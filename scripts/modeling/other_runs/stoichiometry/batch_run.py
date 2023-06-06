import time
import sys
import multiprocessing
import numpy as np

# Runs multiple parallel instances (not communicating with one another) of the given command
# Useful to run multiple parameter-value test runs

SAFETY_N = 10  # Will pause for 1 minutes after 10 processes are started just to stagger the process creation
param_file = sys.argv[1]  # parameters to append to the command as separate arguments
if len(sys.argv) == 3:
    output_file_start = int(sys.argv[2])  # the number to start the output logs at
else:
    output_file_start = 0
param_samples = []  # the actual values of the parameters
param_names = []  # the names of the params
with open(param_file, 'r') as f:
    firstline = True
    first = None
    for line in f:
        if len(line) < 3:
            continue
        if firstline:
            first = line.strip()
            firstline = False
            continue
        line = line.strip()
        param_names.append(line.split(" : ")[0])
        line = line.split(' : ')[1].split(' ')
        param_samples.append(tuple(line))

assert len(np.unique([len(i) for i in param_samples])) == 1, "All parameter samples should be of the same length"
print("Parameters Parsed.")
command_args = first.split(" ")  # the actual command to run


def run_1_instance(args, q, log_file):  # each child process runs this command
    import subprocess
    error_log = log_file.split(".")[0] + "_error.log"
    s = subprocess.run(' '.join(args), stdout=open(log_file, 'w'), stderr=open(error_log, 'w'), shell=True)
    if s.returncode != 0:
        q.put((False, s))
    else:
        q.put((True, s))


if __name__ == '__main__':
    ps = []  # contain all the processes to check after their execution
    q_main = multiprocessing.Queue()  # to return the exit codes of the subprocesses
    for i in range(len(param_samples[0])):
        log_file = "output_log_" + str(i + output_file_start) + ".log"
        param_args = ["-" + param_names[j] + "=" + str(param_samples[j][i]) for j in range(len(param_names))]
        new_command_args = command_args + param_args
        print("BATCH_RUN: RUNNING " + str(new_command_args))
        p = multiprocessing.Process(target=run_1_instance, args=(new_command_args, q_main, log_file), daemon=True)
        p.start()
        if ((i + 1) % SAFETY_N) == 0:
            time.sleep(60)
        ps.append(p)
    print("BATCH_RUN: %d processes have been loaded!" % len(ps))
    for i in range(len(param_samples[0])):  # get all the returned exit codes (wait until all are collected)
        temp = q_main.get()
        print("BATCH_RUN: Recieved 1 return value! ")
        if not temp[0]:
            assert False, "One of the processes returned non-zero error code"
    time.sleep(10)  # to wait for the processes to die (not strictly needed)
    for p in ps:
        if p.is_alive():
            print("Process still alive after all jobs completed. Waiting 10 more seconds")
            time.sleep(10)  # this should not usually do anything more than avert very rare issues
        assert not p.is_alive(), "Process still alive after waiting"
    print("Done!")
