import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import numpy as np
import sys
import os

def rectanlge_1_vs_1(a, b):
    # we want a - b
    clean_boxes = []
    if (a[0] >= b[2])  or (b[0] >= a[2]) or (b[1] >= a[3]) or (a[1] >= b[3]):
        return [a]
    if a[0] < b[0]:
        clean_boxes.append((a[0], a[1], b[0] - 1, a[3]))
    if a[2] > b[2]:
        clean_boxes.append((b[2] + 1, a[1], a[2], a[3]))
    if a[1] < b[1]:
        clean_boxes.append((max(b[0], a[0]), a[1], min(b[2], a[2]), b[1] - 1))
    if a[3] > b[3]:
        clean_boxes.append((max(b[0], a[0]), b[3] + 1, min(b[2], a[2]), a[3]))
    return clean_boxes


def rectangle_all_vs_1(a_list, b):
    clean_boxes = []
    for a in a_list:
        clean_boxes += rectanlge_1_vs_1(a, b)
    return clean_boxes


def rectangle_all_vs_all(a_list, b_list):
    current = a_list[:]
    for b in b_list:
        current = rectangle_all_vs_1(current, b)
    return current

    
def create_patch(x, **kwargs):
    return ptc.Rectangle((x[0], x[1]), x[2] - x[0], x[3] - x[1], **kwargs)


def display_for_example(a_list, b_list, answer=None):
    if answer is None:
        answer = rectangle_all_vs_all(a_list, b_list)
    fig, ax = plt.subplots()
    for a in a_list:
        ax.add_patch(create_patch(a, lw=1, edgecolor='black', facecolor='red', alpha=0.2))
    for b in b_list:
        ax.add_patch(create_patch(b, lw=1, edgecolor='black', facecolor='none', hatch='/'))
    for c in answer:
        ax.add_patch(create_patch(c, lw=1, edgecolor='black', facecolor='none', hatch='.'))
    plt.xlim((np.min(np.array([x[0] for x in a_list + b_list]).flatten()) - 1, np.max(np.array([x[2] for x in a_list + b_list]).flatten()) + 1))
    plt.ylim((np.min(np.array([x[1] for x in a_list + b_list]).flatten()) - 1, np.max(np.array([x[3] for x in a_list + b_list]).flatten()) + 1))


def create_boxes_from_ranges(x, y):
    x = x.split('),')
    y = y.split('),')
    x = [i.replace('(', '').replace(')', '').replace(' ', '').split(',') for i in x]
    y = [i.replace('(', '').replace(')', '').replace(' ', '').split(',') for i in y]
    clean_boxes = []
    for x_sub in x:
        for y_sub in y:
            clean_boxes.append((int(x_sub[0]), int(y_sub[0]), int(x_sub[1]), int(y_sub[1])))
    return clean_boxes

def reverse_of_above(x):
    return f'{x[0]} - {x[2]} : {x[1]} - {x[3]}'

#display_for_example([(0, 0, 1 , 10)], [(2, 1, 9, 3)])
#display_for_example([(0, 0, 10 , 10)], [(11, 11, 12, 12)])
#display_for_example([(0, 0, 10 , 10)], [(11, 11, 12, 12), (2, 1, 9, 3)])
#display_for_example([(0, 0, 10 , 10)], [(1, 1, 2, 9), (11, 11, 12, 12)])
#display_for_example([(0, 0, 10 , 10)], [(1, 1, 2, 9), (2, 1, 9, 3)])

path = sys.argv[1]
names_cutoff = [(x.split('_')[2], x.split('_')[3], int(x.split('_')[4].replace('.txt', ''))) for x in os.listdir(path) if 'list_contacts' in x]
unique_names = list(set(([x[:2] for x in names_cutoff])))   
for n1, n2 in unique_names:
    f1 = os.path.join(path, f'list_contacts_{n1}_{n2}_20.txt')
    f2 = os.path.join(path, f'list_contacts_{n1}_{n2}_25.txt')
    with open(f1) as f:
        rd = f.read()
    rd = rd.strip().split('\n')[1:]
    rd = [i.split('\t') for i in rd]
    a_list = []
    for i in rd:
        a_list += create_boxes_from_ranges(i[0], i[1])
    with open(f2) as f:
        rd = f.read()
    rd = rd.strip().split('\n')[1:]
    rd = [i.split('\t') for i in rd]
    b_list = []
    for i in rd:
        b_list += create_boxes_from_ranges(i[0], i[1])
    if (len(a_list) == 0) or (len(b_list) == 0):
        continue
    if (n1 == 'GPKP1a') or (n2 == 'GPKP1a'):
        continue
    answer = rectangle_all_vs_all(a_list, b_list)
    display_for_example(a_list, b_list, answer)
    plt.title(f'{n1} -> {n2}')
    plt.show()
    print(f'{n1} -> {n2}')
    for ans in answer:
        print('\t', reverse_of_above(ans))
