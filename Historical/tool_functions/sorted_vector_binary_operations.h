#pragma once

template <typename T>
int sorted_vector_binary_operations(std::vector<T> &input_vector, int key)
{

	/*return the possible insertion position*/

	int left = 0, right = input_vector.size() - 1;
	int result = -1;
	while (left <= right) // it will be skept when input_vector.size() == 0
	{
		int mid = (left + right) / 2; // mid is between left and right (may be equal);
		if (input_vector[mid] == key)
		{
			result = mid;
			left = mid + 1;
		}
		if (input_vector[mid] > key)
		{
			right = mid - 1; // the elements after right are always either empty, or have larger keys than input key
		}
		else
		{
			left = mid + 1; // the elements before left are always either empty, or have smaller keys than input key
		}
	}
	return result == -1 ? left : result;
};