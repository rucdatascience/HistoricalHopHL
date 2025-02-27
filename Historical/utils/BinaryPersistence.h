#pragma once

#include <fstream>
#include <vector>
namespace experiment
{

	template <typename T, typename Enable = void>
	class BinarySerializer
	{
	public:
		static void saveBinary(std::ofstream &out, const T &data)
		{
			out.write(reinterpret_cast<const char *>(&data), sizeof(T));
		}

		static void loadBinary(std::ifstream &in, T &data)
		{
			in.read(reinterpret_cast<char *>(&data), sizeof(T));
		}
	};

	template <typename K, typename V>
	class BinarySerializer<std::pair<K, V>>
	{
	public:
		static void saveBinary(std::ofstream &out, const std::pair<K, V> &pair)
		{
			BinarySerializer<K>::saveBinary(out, pair.first);
			BinarySerializer<V>::saveBinary(out, pair.second);
		}

		static void loadBinary(std::ifstream &in, std::pair<K, V> &pair)
		{
			BinarySerializer<K>::loadBinary(in, pair.first);
			BinarySerializer<V>::loadBinary(in, pair.second);
		}
	};

	template <typename T>
	class BinarySerializer<std::vector<T>>
	{
	public:
		static void saveBinary(std::ofstream &out, const std::vector<T> &vec)
		{
			size_t size = vec.size();
			BinarySerializer<size_t>::saveBinary(out, size);
			for (const auto &item : vec)
			{
				BinarySerializer<T>::saveBinary(out, item);
			}
		}

		static void loadBinary(std::ifstream &in, std::vector<T> &vec)
		{
			size_t size;
			BinarySerializer<size_t>::loadBinary(in, size);
			vec.resize(size);
			for (auto &item : vec)
			{
				BinarySerializer<T>::loadBinary(in, item);
			}
		}
	};

	template <typename T>
	inline void saveBinary(std::ofstream &out, const T &data)
	{
		BinarySerializer<T>::saveBinary(out, data);
	}

	template <typename T>
	inline void loadBinary(std::ifstream &in, T &data)
	{
		BinarySerializer<T>::loadBinary(in, data);
	}
}