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
# work_2 = [ "","","","" ]
work_3 = {"item1":1, "item2":3}
# work_3 = {"item2":3}

def multiproc_cmd_list_tuple(arg_list1, arg_list2, **kwargs):
    '''working doc-string'''

    # Empty list
    c_list = list()

    # for idx,file in enumerate(arg_list1):
    for idx in range(0, len(arg_list1), 1):
        tmp_list = [arg_list1[idx], arg_list2[idx]]
        # print(tmp_list)
        for key, item in kwargs.items():
            # tmp_dict = {key: item}
            # new_dict.update(tmp_dict)
            # tmp_list.append(f"{key}={item}")
            tmp_list.append(f"{item}")
        c_list.append(tmp_list)

    c_tup = tuple(c_list)
    return c_tup

def print_func(n,a,item1="SomeItem",item2="MoreItems"):
    print(n,a,item1,item2)

if __name__ == "__main__":
    cmd_tup = multiproc_cmd_list_tuple(work_1,work_2,**work_3)
    # cmd_tup = multiproc_cmd_list_tuple(work_1, work_2)
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.starmap(print_func, cmd_tup)
    print(results)