#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class HashTable {
	public:
		static int STRING_HASH_CONSTANT;
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
		void* remove(int intKey);
		bool insert(void* key, int keySize, void* entry);
		void* find(void* key, int keySize);
		void* remove(void* key, int keySize);
	
	private:
		int baseHash(int size, int hashKey);
		int stepHash(int size, int hashKey);
		void doubleSize();
		int getIntKey(char* inputString, int len);


};
