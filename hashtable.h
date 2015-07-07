#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class HashTable {
	public:
		static double BASE_HASH_CONSTANT;
		static double STEP_HASH_CONSTANT;
		static int INITIAL_TABLE_SIZE;

		int* dummy;
		int* keys;
		void** entries;
		int size,load;

		HashTable();
		~HashTable();
		bool insert(int intKey,void* entry);
		void* find(int intKey);
		bool remove(int intKey);
	
	private:
		int baseHash(int size, int hashKey);
		int stepHash(int size, int hashKey);
		void doubleSize();


};
