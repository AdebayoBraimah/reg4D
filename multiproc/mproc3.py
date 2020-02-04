import multiprocessing
import time

work_1 = [ "A", "B", "C", "D" ]
work_2 = [ 5, 2, 1, 3 ]

# def work_log(work_1,work_2):
#     print(f"process {work_1}, {work_2}")
#     time.sleep(int(work_2))
#     print(f"Process {work_1} has finished")
#
# def pool_handler():
#     p = multiprocessing.Pool(2)
#     p.map(work_log(work_1,work_2),work_2)
#
# if __name__ == "__main__":
#     pool_handler()

work_1 = [ "A", "B", "C", "D" ]
work_2 = [ 5, 2, 1, 3 ]

def print_func(n,a):
    print(n,a)

if __name__ == "__main__":
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.starmap(print_func, (work_1,work_2))
    print(results)