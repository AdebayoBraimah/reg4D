import multiprocessing
import time
import sys

# multiproc example
def print_func(continent="Asia"):
    '''Test func for multi-proc'''
    # print(f"The name of continent is: {continent}")
    print("The name of continent is: ", continent)


if __name__ == "__main__":
    names = ['America', 'Europe', 'Africa']
    procs = []
    # proc = multiprocessing.Process(target=print_func)

    for name in names:
        # print(name)
        proc = multiprocessing.Process(target=print_func, args=(name,))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()