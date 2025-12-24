import math

# recursive version
def partitions(n, k):
    # Optimized further by adding k == 1 as a base case because there is always only 1 way to partition with 1 part
    if n == 0 or k == 1:
        return 1
    if n < 0:
        return 0
    return partitions(n, k - 1) + partitions(n - k, k)

print(partitions(22, 6)) 

# can we optimize furter/not use recursion? Not really
# def partitions_new(n, k):
#     sum = 0
#     i = 0
#     for i in range(n, 0, -k):
#         sum += math.floor(i/k)+1
#         print("sum" + str(sum))
#     return sum

# print(partitions_new(10, 3))


