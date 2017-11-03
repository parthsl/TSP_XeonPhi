print("Enter filename: ")
filename = raw_input()
print("Enter number of cities: ")
city = int(raw_input())

city = range(city)

with open(filename) as f:
    array = [[int(x) for x in line.split()] for line in f][0]

print("Unecessary cities:")
print( list(set(array)-set(city)))
print("Not in list:")
print( list(set(city)-set(array)))
