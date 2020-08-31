# import multiprocessing
# import time
# import sys
print('hello world')


# start = time.perf_counter()


# def do_something():
#     print('sleeping 1 second(s)...')
#     time.sleep(1)
#     print('Done sleeping')
#
#
# processes = []  # create a list to store these submitted processes
# for _ in range(0, 10):
#     p = multiprocessing.Process(target=do_something)
#     p.start()
#     processes.append(p)
# for process in processes:
#     process.join()
#
# finish = time.perf_counter()  # submit the next process, which is the same as the synchronous process
# print('Finished in {} sec'.format(finish - start))