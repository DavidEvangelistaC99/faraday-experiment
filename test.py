import numpy
a = [[1, 2, 4], [2, 3, 3]]

arr = numpy.array(a, dtype='int')
print(arr.shape)

for i in range(int(len(arr[0]))):
    print(arr[:,i])

xx = numpy.zeros([4, 4])

print(xx)