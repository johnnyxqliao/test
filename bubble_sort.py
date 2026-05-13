import time


def bubble_sort(numbers):
    """Use bubble sort to return numbers in ascending order."""
    result = numbers[:]
    length = len(result)

    for i in range(length):
        swapped = False
        for j in range(0, length - i - 1):
            if result[j] > result[j + 1]:
                result[j], result[j + 1] = result[j + 1], result[j]
                swapped = True

        if not swapped:
            break

    return result


if __name__ == "__main__":
    nums = [37, 82, 14, 96, 5, 63, 28, 71, 44, 19, 90, 7, 56, 33, 68, 12, 100, 41]

    start_time = time.perf_counter()
    sorted_nums = bubble_sort(nums)
    end_time = time.perf_counter()

    print("Before sort:", nums)
    print("After sort:", sorted_nums)
    print(f"Run time: {end_time - start_time:.8f} seconds")
