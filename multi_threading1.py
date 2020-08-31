import concurrent.futures
import time

start = time.perf_counter()


def do_something(seconds):
    print(f'Sleeping for {seconds} second(s)')
    time.sleep(seconds)
    return f'Done sleeping for {seconds} second(s)'


secs = [5, 4, 3, 2, 1]

with concurrent.futures.ThreadPoolExecutor() as executor:
    threads = [executor.submit(do_something, sec) for sec in secs]  # map will return results list, submit will return a future object
    for f in concurrent.futures.as_completed(threads):
        print(f.result())

finish = time.perf_counter()
print(f'Finished in {finish - start} sec(s)')
