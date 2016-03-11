#ifndef CM_HEAP_WRAPPER_H
#define CM_HEAP_WRAPPER_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

//DESCRIPTION:
//Template class for allocation of (numerous) small objects
//Allocates objects in blocks to avoid memory overhead
//Allocated memory can be reused -> reset()


// new HeapWrapper<_simple_rkmer>( (_4mb/sizeof(_simple_rkmer)) -1 );
template <class T> class HeapWrapper{

	public:
		typedef T Object_type;
		class _block{
			public:
				_block      *next; 
				size_t      index;
				Object_type *data; 
				_block(	const size_t blocksize, 
							const size_t index, 
							_block *next=0):
							next(next), 
							index(index){ 
					data = new Object_type[blocksize];
				}
				~_block(){ delete [] data; }	
		};

		HeapWrapper( const size_t blocksize ):_BLOCK_SIZE(blocksize){
			_first   = new _block( _BLOCK_SIZE, 0 );
			_current = _first;
			_i       = 0;
			_bc      = 1;
		}

		virtual ~HeapWrapper(){
			while(_first!=0){
				_current = _first;
				_first = _first->next;
				delete _current;
			}
		}

		void reset(){
			_current = _first;
			_i       = 0;
/*			while(_first!=0){
				_current = _first;
				_first = _first->next;
				delete _current;
			}
			_first   = new _block( _BLOCK_SIZE, 0 );
			_current = _first;
			_i       = 0;
			_bc      = 1;*/

		}

		inline Object_type* get_next(){
			if( _i==_BLOCK_SIZE ){
				_block *tmp = _current->next;
				if( tmp==0){
					tmp            = new _block(_BLOCK_SIZE, _bc++);
					_current->next = tmp;
				}
				_current = tmp;
				_i       = 0;
			}
			return &_current->data[_i++];
		}

		inline Object_type* get_first(){ return &_first->data[0]; }

		size_t get_block_count(){ return _bc; }

		size_t get_current_memory_usage(){ 
			return (_current->index*_BLOCK_SIZE*sizeof(Object_type)) +  (_i*sizeof(Object_type)); 
		}

		size_t get_number_of_current_objects_in_use(){ return _current->index*_BLOCK_SIZE + _i; }

		size_t get_instance_memory_usage(){
			return ((_bc*_BLOCK_SIZE*sizeof(Object_type)) + 4*_bc*sizeof(size_t) + _bc*sizeof(_block));
		}


	protected:
		_block *_first, *_current;
		const size_t _BLOCK_SIZE;
		size_t _i;
		size_t _bc;
};

#endif
