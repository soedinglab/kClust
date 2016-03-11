#ifndef CM_SIMPLE_LIST_H
#define CM_SIMPLE_LIST_H
/***************************************************************************
 *   Copyright (C) 2006 by Christian Mayer                                 *
 *   christian.eberhard.mayer@googlemail.com                               *
 ***************************************************************************/

//DESCRIPTION:
//Template class implementing a single linked list
//std::list is double linked

#include <cstdlib>
template <class T> class Simple_list{

	public:
		typedef T Element_value_type;
			
		class Element{
			public:
				Element();
				Element(const Element_value_type&, Element *);
				Element(const Element&);
				~Element();
				Element *next;
				Element_value_type value;
		};

		Simple_list();
		Simple_list(const Element_value_type&);
		~Simple_list();

		void append(const Element_value_type&);
			
		class Iterator{
			public:
				Iterator(Element *e=0):e(e){}
				~Iterator(){}
				Iterator& operator=(const Iterator &it){
					e = it.e;
					return (*this);
				}
				bool operator==(const Iterator &it){
					return(it.e==e);
				}
				bool operator!=(const Iterator& it){
					return(it.e!=e);
				}
				Iterator& operator++(){
					e=e->next;
					return (*this);
				}
				Iterator operator++(int){
					Iterator tmp(*this);
					++(*this);
					return tmp;
				}
				Element_value_type& operator*(){
					return e->value;	
				}
			private:
				Element *e;
		};


		Iterator begin(){ return _first->next; }
		Iterator end()  { return 0;            }

		inline size_t length() const { return len; }

	private:
		Element *_first, *_last;
		size_t len;
};


template <class T> Simple_list<T>::Simple_list(){
	_first  = new Element();
	_last = _first;
	len=0;
}

template <class T> Simple_list<T>::Simple_list(const Element_value_type &value){
	_first  = new Element();
	_last = _first;
	len=0;
	append(value);
}

template <class T> Simple_list<T>::~Simple_list(){
	while(_first!=0){
		_last = _first;
		_first = _first->next;
		delete _last;	
	}
}

template <class T> void Simple_list<T>::append(const T &value){
	_last->next = new Element(value);
	_last       = _last->next;
	++len;
}

template <class T> Simple_list<T>::Element::Element():next(0){}

template <class T> Simple_list<T>::Element::Element(const Element_value_type &value, Element *next=0):next(next),value(value){}

template <class T> Simple_list<T>::Element::Element(const Element &e){
	next = e->next;
	value( e.value );
}

template <class T> Simple_list<T>::Element::~Element(){}

#endif
