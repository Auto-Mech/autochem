from autoreact import params

test_arr = [1e12, 1.5, 50000]
test_arr2 = [1e12, 1.5, 50000]

params = params.RxnParams()
params.set_arr(test_arr, test_arr2)

print(params.arr)
