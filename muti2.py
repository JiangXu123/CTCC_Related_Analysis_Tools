import concurrent.futures
import time


start = time.perf_counter()


def do_something(seconds):
    print(f'sleeping {seconds} second(s)...')
    time.sleep(seconds)
    return f'Done sleeping in {seconds} seconds(s)'
secs = [5 ,4, 3, 2, 1]

with concurrent.futures.ProcessPoolExecutor() as executor: # do the same thing without list comprehension
    results = [executor.submit(do_something, sec) for sec in secs]
    # The above code return a list of future object, stored in results
    for f in concurrent.futures.as_completed(results):  # using the list comprehension, it's possible to use the as_completed method
        print(f.result())

finish = time.perf_counter()    # submit the next process, which is the same as the synchronous process
print('Finished in {} sec'.format(finish-start))