from IMP.isd import create_gmm
import time
import copy
import multiprocessing


# Generate multiple EM maps using Gaussian Mixtures with differing number of centers to find the smallest
# mixture satisfying a cross correlation cutoff

def worker_foo(q):
    from IMP.isd import create_gmm
    while not q.empty():  # weaker check. May Block the next instruction.
        arg_obj = q.get()
        create_gmm.run(arg_obj)
        q.task_done()


if __name__ == '__main__':
    # read the arguments for the create_gmm script
    range_size = 60  # range of the n_centers for GMM tried
    args = create_gmm.parse_args()
    start_time = time.time()

    # alter the n-centers, output file names
    main_out_file = args.out_file
    main_out_map = args.out_map
    start_val = args.n_centers
    args_list = []
    for i in range(start_val, start_val + range_size, 10):
        args.n_centers = i
        args.out_file = main_out_file + str(i) + '.txt'
        args.out_map = main_out_map + str(i) + '.mrc'
        args_list.append(copy.deepcopy(args))

    q_main = multiprocessing.JoinableQueue()
    n_processes = 6
    for i in args_list:
        q_main.put(i)
    for i in range(n_processes):
        multiprocessing.Process(target=worker_foo, daemon=True, args=(q_main,)).start()
    q_main.join()  # wait for the tasks to be done
    print('Total time taken: ' + str(round(time.time() - start_time, 2)) + ' seconds.')
    # The stdout will be a jumble due to multi-process writing to STDOUT

    # Cross-correlation of all of the generated .mrc with the original density using Chimera
    # Command file: compare_gmm.com and the output stored as output_compare_gmm.txt
