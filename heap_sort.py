import time


def heapify(numbers, length, root_index):
    largest = root_index
    left_child = 2 * root_index + 1
    right_child = 2 * root_index + 2

    if left_child < length and numbers[left_child] > numbers[largest]:
        largest = left_child

    if right_child < length and numbers[right_child] > numbers[largest]:
        largest = right_child

    if largest != root_index:
        numbers[root_index], numbers[largest] = numbers[largest], numbers[root_index]
        heapify(numbers, length, largest)


def heap_sort(numbers):
    """Use heap sort to return numbers in ascending order."""
    result = numbers[:]
    length = len(result)

    for i in range(length // 2 - 1, -1, -1):
        heapify(result, length, i)

    for i in range(length - 1, 0, -1):
        result[0], result[i] = result[i], result[0]
        heapify(result, i, 0)

    return result


if __name__ == "__main__":
    nums = [37, 82, 14, 96, 5, 63, 28, 71, 44, 19, 90, 7, 56, 33, 68, 12, 100, 41, 25, 79, 3, 88, 52, 16, 74]

    start_time = time.perf_counter()
    sorted_nums = heap_sort(nums)
    end_time = time.perf_counter()

    print("Before sort:", nums)
    print("After sort:", sorted_nums)
    print(f"Run time: {end_time - start_time:.8f} seconds")
