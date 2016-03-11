#ifndef CM_SIMPLE_HASH_H
#define CM_SIMPLE_HASH_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//Template class implementing a very basic hash

#include <cstddef>
#include <cmath>
#include <inttypes.h>
#include <string.h>

#define MIX(a,b,c){\
  a=a-b;  a=a-c;  a=a^(c>>13); \
  b=b-c;  b=b-a;  b=b^(a<<8);  \
  c=c-a;  c=c-b;  c=c^(b>>13); \
  a=a-b;  a=a-c;  a=a^(c>>12); \
  b=b-c;  b=b-a;  b=b^(a<<16); \
  c=c-a;  c=c-b;  c=c^(b>>5);  \
  a=a-b;  a=a-c;  a=a^(c>>3);  \
  b=b-c;  b=b-a;  b=b^(a<<10); \
  c=c-a;  c=c-b;  c=c^(b>>15); \
}
// 0.618. * 2^32
#define GOLDEN_RATIO 0x9e3779b9

template <class T> class Simple_hash{

	public:
		typedef uint32_t       u4;   /* unsigned 4-byte type */
		typedef unsigned char  u1;   /* unsigned 1-byte type */
		typedef T value_type;
	
		class Hash_entry{
			public:
				Hash_entry():	key(0), 
									value(0), 
									hash(0), 
									next(0){}
				Hash_entry(	const char * const k,
								u1 len, 
								u4 h, 
								value_type v, 
								Hash_entry *last, 
								Hash_entry* n=0){
									key = new char[len+1];
									memcpy(key, k, sizeof(char)*len);
									key[len] = '\0';
									hash  = h;
									value = v;
									next  = n;
									last_ins = last;
				}
				~Hash_entry(){ delete [] key; }	
				char       *key;
				u4         hash;
				value_type value;
				Hash_entry *next;
				Hash_entry *last_ins;
		};

		class Iterator{
			public:
				friend class Simple_hash<value_type>;
				Iterator(Hash_entry* entry=0):entry(entry){}
				~Iterator() {}
				Iterator& operator=(const Iterator& other){
					entry = other.entry;
					return (*this);
				}
				bool operator==(const Iterator& other){ return (entry == other.entry); }
				bool operator!=(const Iterator& other){ return (entry != other.entry); }
				Iterator& operator++(){
					if( entry!=0 )  entry = entry->last_ins;	
					return(*this);
				}
				Iterator operator++(int){
					Iterator tmp(*this);
					++(*this);
					return(tmp);
				}
				char* operator*()      { return entry->key; }
				value_type& getValue() { return entry->value;}
			private:
				Hash_entry *entry;
		};

		Simple_hash(u4 size){
			size_t bits = int( log10(size)/log10(2.0f) )+1;
			data_size   = 1 << bits;
			mask        = data_size-1;
			data        = new Hash_entry*[data_size];
			memset( data, 0, sizeof(Hash_entry*)*data_size );	
			last_ins    = 0;
			mem         = 8 + data_size*sizeof(value_type*);
			load        = 0;
		}

		~Simple_hash(){
			Hash_entry *tmp;
			while(last_ins!=0){
				tmp  = last_ins;
				last_ins = last_ins->last_ins;
				delete tmp;
			}
			delete [] data;	
		}
		
		// key+len = KEY, value = VALUE
		void put( const char * const key, const u1 len, value_type value){
			u4 h = _hash((const u1*)key, len);
			Hash_entry *tmp = data[h];
			if( tmp==0 ){
				last_ins = new Hash_entry(key, len, h, value, last_ins);
				data[h]  = last_ins;
				++load;
				return;
			}
			while( tmp!=0 ){ 
				bool b=true;	
				for(u1 i=0; i<len; ++i){
					if( tmp->key[i]!=key[i] ){
						b=false;
						break;
					}
				}
				if(b){
					tmp->value = value;
					return;
				}
				if( tmp->next==0 ){
					last_ins  = new Hash_entry(key, len, h, value, last_ins);
					tmp->next = last_ins;
					++load;
					break;
				}
				tmp  = tmp->next;
			}
		}
		
		// key + len = KEY
		value_type* get(const char * const key, u1 len){
			u4 h = _hash((const u1*)key, len);
			Hash_entry *tmp  = data[h];			
			while( tmp!=0 ){ 
				bool b=true;	
				for(u1 i=0; i<len; ++i){
					if( tmp->key[i]!=key[i] ){
						b=false;
						break;
					}
				}
				if(b) return &(tmp->value);
				tmp  = tmp->next;
			}
			return 0;
		}

		bool remove(const char * const key, u1 len){
			u4 h = _hash((const u1*)key, len);
			Hash_entry *tmp  = data[h];
			
			if( tmp==0 ) return false;
			Hash_entry *prev_in_bucket = 0;
			while( tmp!=0 ){ 
				bool b=true;	
				for(u1 i=0; i<len; ++i){
					if( tmp->key[i]!=key[i] ){
						b=false;
						break;
					}
				}
				if(b){
					if( tmp==last_ins ) last_ins = tmp->last_ins;
					if( prev_in_bucket==0 ){
						data[h] = tmp->next;
						delete tmp; 
					}else{
						prev_in_bucket->next = tmp->next;
						delete tmp;
					}
					--load;
					return true;
				}
				prev_in_bucket = tmp;
				tmp  = tmp->next;
			}
			return false;
		}

		bool exists(const char * const key, u1 len){
			u4 h = _hash((const u1*)key, len);
			Hash_entry *tmp  = data[h];			
			if( tmp==0 ) return false;
			while( tmp!=0 ){ 
				bool b=true;	
				//printf("CHECK:\n");
				for(u1 i=0; i<len; ++i){
					//printf("     %i %i \n",tmp->key[i],key[i]);
					if( tmp->key[i]!=key[i] ){
						b=false;
						break;
					}
				}
				if(b) return true;
				tmp  = tmp->next;
				//printf("TMP IS NOW: %x\n", tmp);
			}
			return false;
		}

		Iterator begin(){
			return Iterator(last_ins);
		}

		Iterator end(){
			return Iterator(0);
		}

		inline const size_t get_load() const{ return load; }

		const size_t get_memory_usage() const{ 
			return mem + ( (sizeof(Hash_entry)+sizeof(Hash_entry*)+sizeof(size_t)) * load ); 
		}

	private:
		u4 data_size;
		u4 load;
		u4 mask;
		size_t mem;

		Hash_entry **data;
		Hash_entry *last_ins;

		//very dispersive hash function of Bob Jenkins from http://www.burtleburtle.net/bob/hash/doobs.html
		u4 _hash(const u1 *k, u1 length){
   			register u4 a,b,c;     /* the internal state */
   			u4          len;       /* how many key bytes still need mixing */
   			/* Set up the internal state */
   			len = length;
   			a = b = GOLDEN_RATIO;  /* the golden ratio; an arbitrary value */
   			c = 0xabcd;           /* variable initialization of internal state */
   			/*---------------------------------------- handle most of the key */
			while (len >= 12){
			a=a+(k[0]+((u4)k[1]<<8)+((u4)k[2]<<16) +((u4)k[3]<<24));
			b=b+(k[4]+((u4)k[5]<<8)+((u4)k[6]<<16) +((u4)k[7]<<24));
			c=c+(k[8]+((u4)k[9]<<8)+((u4)k[10]<<16)+((u4)k[11]<<24));
			MIX(a,b,c);
			k = k+12; len = len-12;
			}			
			/*------------------------------------- handle the last 11 bytes */
			c = c+length;
			switch(len){
				case 11: c=c+((u4)k[10]<<24);
				case 10: c=c+((u4)k[9]<<16);
				case 9 : c=c+((u4)k[8]<<8);
				/* the first byte of c is reserved for the length */
				case 8 : b=b+((u4)k[7]<<24);
				case 7 : b=b+((u4)k[6]<<16);
				case 6 : b=b+((u4)k[5]<<8);
				case 5 : b=b+k[4];
				case 4 : a=a+((u4)k[3]<<24);
				case 3 : a=a+((u4)k[2]<<16);
				case 2 : a=a+((u4)k[1]<<8);
				case 1 : a=a+k[0];
				/* case 0: nothing left to add */
			}
			MIX(a,b,c);
			return c&mask;
		}
};

#endif
