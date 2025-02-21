#pragma once
// 通用二进制序列化模板
#include <iostream>
#include <fstream>
#include <type_traits>
#include <vector>
namespace experiment {
	template <typename T>
	typename std::enable_if<std::is_trivially_copyable<T>::value, void>::type
		saveBinary(std::ofstream& out, const T& data) {
		out.write(reinterpret_cast<const char*>(&data), sizeof(T));
	}

	template <typename T>
	typename std::enable_if<std::is_trivially_copyable<T>::value, void>::type
		loadBinary(std::ifstream& in, T& data) {
		in.read(reinterpret_cast<char*>(&data), sizeof(T));
	}

	template <typename K, typename V>
	void saveBinary(std::ofstream& out, const std::pair<K, V>& pair) {
		saveBinary(out, pair.first);
		saveBinary(out, pair.second);
	}

	template <typename K, typename V>
	void loadBinary(std::ifstream& in, std::pair<K, V>& pair) {
		loadBinary(in, pair.first);
		loadBinary(in, pair.second);
	}

	template <typename T>
	void saveBinary(std::ofstream& out, const std::vector<T>& vec) {
		size_t size = vec.size();
		saveBinary(out, size);
		for (const auto& item : vec) {
			saveBinary(out, item);
		}
	}

	template <typename T>
	void loadBinary(std::ifstream& in, std::vector<T>& vec) {
		size_t size;
		loadBinary(in, size);
		std::vector<T>().swap(vec);
		vec.resize(size);
		for (auto& item : vec) {
			loadBinary(in, item);
		}
	}
}
