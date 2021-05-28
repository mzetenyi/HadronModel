#ifndef MULTIARRAY_HPP
#define MULTIARRAY_HPP

#include <vector>
#include <map>
#include <iostream>
//#include <cstdarg>
#include <cassert>
#include <initializer_list>
#include "halfint.hpp"
#include "utils.hpp"

using namespace std;

class idx {
public:
  idx();
  idx(int);
  idx(halfint _min, halfint _max, halfint _delta=_1);
  halfint min() const;
  halfint max() const;
  halfint delta() const;
  int dim() const;
  int operator()(halfint h) const;
  bool match(halfint h) const;
  friend bool operator==(const idx& lhs, const idx& rhs);
  friend ostream& operator<<(ostream&, const idx&);
  static const idx idx_null;
  static const idx idx_0;
  static const idx idx_lor;
  static const idx idx_s1;
  static const idx idx_s10;
  static const idx idx_s1h;
  static const idx idx_s3h;
  static const idx idx_s5h;
private:
  halfint _min;
  halfint _max;
  halfint _delta;
  int _dim;
};

bool operator!=(const idx& lhs, const idx& rhs);

extern const idx& idx_null;  // placeholder for non-existing index
extern const idx& idx_0;     // one element array  
extern const idx& idx_lor;   // lorentz index
extern const idx& idx_s1;    // spin 1             
extern const idx& idx_s1_0;  // spin 1 massless    
extern const idx& idx_s1h;   // spin 1/2           
extern const idx& idx_s3h;   // spin 3/2           
extern const idx& idx_s5h;   // spin 5/2           

const idx& idx_spin(halfint spin);

class NamedIndex {
public:
  NamedIndex() = delete;
  NamedIndex(const vector<string>& names);
  halfint operator()(const string& name) const;
  string operator()(int i) const;
  vector<string> getNames() const;
  idx getIdx() const;
  static vector<idx> getIdxList(const vector<NamedIndex> NIs);
  friend ostream& operator<<(ostream&, const NamedIndex&);
private:
  vector<string> names;
  map<string,int> name_to_int;
  map<int,string> int_to_name;
};

/**
   Represent an N dimensional array of elements of type T.

   \todo The dimension and the index types of the array is not known at declaration only at 
   construction. i.e. at runtime. Therefore there's no compile time check on the indices
   of the array.
*/
template <typename T>
class MultiArray {
public:
  /**
     Constructor. The array is filled with default constructed
     elements. This depends on the implementation of the vector.resize()
     method, therefore there is no initialization for primitive types
     like int or double.
  */
  MultiArray(vector<idx> indices) :
    _idx(indices) {
    nelement = 1;
    for (idx i : _idx) {
      nelement *= i.dim();
    }
    vals.resize(nelement);
  }

  MultiArray(vector<NamedIndex> IDs) :
    _idx(NamedIndex::getIdxList(IDs)) {
    nelement = 1;
    for (idx i : _idx) {
      nelement *= i.dim();
    }
    vals.resize(nelement);
  }

  MultiArray(idx _idx0, idx _idx1=idx_null, idx _idx2=idx_null, idx _idx3=idx_null, idx _idx4=idx_null,
             idx _idx5=idx_null, idx _idx6=idx_null, idx _idx7=idx_null, idx _idx8=idx_null, idx _idx9=idx_null) {
    vector<idx> args{_idx0, _idx1, _idx2, _idx3, _idx4, _idx5, _idx6, _idx7, _idx8, _idx9};
    for (auto it=args.begin(); *it != idx_null; it++) {
      _idx.push_back(*it);
    }
    nelement = 1;
    for (idx i : _idx) {
      nelement *= i.dim();
    }
    vals.resize(nelement);
  }

  void fill(const T& val) {
    for (int i(0); i<nelement; i++) { vals[i] = val; }
  }

  T& operator()(vector<halfint> indices) {
    //    cerr << "in T& MultiArray::operator()(vector<halfint> indices)" << endl;
    // cerr << "called for :" << endl;
    // cerr << *this << endl;
    // cerr << "---------------------------------" << endl;
    // cerr << "indices: ";
    // for (auto index : indices) cerr << index << " ";
    // cerr << endl;
    vector<int> i(indices.size());
    for (int k(0); k<indices.size(); k++) {
      i[k] = _idx[k](indices[k]);
      //      cerr << "k, i[k]: " << k << " " << i[k] << endl;
    }
    int ind = index(i);
    // cerr << "ind: " << ind << endl;
    // cerr << "nelement: " << nelement << endl;
    //    cerr << "return from T& MultiArray::operator()(vector<halfint> indices)" << endl;
    return vals[ind];
  }

  const T& operator()(vector<halfint> indices) const {
    //    cerr << "in const T& MultiArray::operator()(vector<halfint> indices) const" << endl;
    vector<int> i(indices.size());
    for (int k(0); k<indices.size(); k++) {
      i[k] = _idx[k](indices[k]);
    }
    int ind = index(i);
    //    cerr << "return from const T& MultiArray::operator()(vector<halfint> indices) const" << endl;
    return vals[ind];
  }

  T& operator()(halfint h0, halfint h1=_0, halfint h2=_0, halfint h3=_0, halfint h4=_0,
                halfint h5=_0, halfint h6=_0, halfint h7=_0, halfint h8=_0, halfint h9=_0) {
    vector<halfint> args{h0, h1, h2, h3, h4, h5, h6, h7, h8, h9};
    args.resize(_idx.size());
    return operator()(args);
  }
  
  const T& operator()(halfint h0, halfint h1=_0, halfint h2=_0, halfint h3=_0, halfint h4=_0,
                halfint h5=_0, halfint h6=_0, halfint h7=_0, halfint h8=_0, halfint h9=_0) const {
    vector<halfint> args{h0, h1, h2, h3, h4, h5, h6, h7, h8, h9};
    args.resize(_idx.size());
    return operator()(args);
  }

  T& operator()(int i0, int i1=0, int i2=0, int i3=0, int i4=0,
                int i5=0, int i6=0, int i7=0, int i8=0, int i9=0) {
    vector<halfint> args{halfint(i0), halfint(i1), halfint(i2), halfint(i3), halfint(i4),
        halfint(i5), halfint(i6), halfint(i7), halfint(i8), halfint(i9)};
    args.resize(_idx.size());
    return operator()(args);
  }

  const T& operator()(int i0, int i1=0, int i2=0, int i3=0, int i4=0,
                int i5=0, int i6=0, int i7=0, int i8=0, int i9=0) const {
    vector<halfint> args{halfint(i0), halfint(i1), halfint(i2), halfint(i3), halfint(i4),
        halfint(i5), halfint(i6), halfint(i7), halfint(i8), halfint(i9)};
    args.resize(_idx.size());
    return operator()(args);
  }

  vector<idx> get_idx() const {
    return _idx;
  }

  friend ostream& operator<<(ostream& out, const MultiArray& arr) {
    out << "MultiArray object with indices:" << endl;
    for (auto index : arr._idx) cerr << index << endl;
    return out;
  }
  
private:
  vector<idx> _idx;
  int nelement;
  vector<T> vals;

  // default construction forbidden:
  MultiArray();

  // copy construction forbidden:
  //  Array(const Array<T,N>&);

  int index(const vector<int>& i) const {
    /*
    cerr << "Array::index" << endl;
    cerr << "( ";
    for (int k(0); k<N; k++) { cerr << i[k] << ", "; }
    cerr << ")   of   ( ";
    for (int k(0); k<N; k++) { cerr << _idx[k].dim() << ", "; }
    cerr << ")" << endl;
    //*/
    for (int k(0); k<i.size(); k++) {
      if (!(i[k]<_idx[k].dim())) {
        cerr << "something wrong in MultiArray object" << endl;
        cerr << *this << endl;
        PR(i.size());
        for (int ii(0); ii<i.size(); ii++) {
          PL(ii); PL(i[ii]); PR(_idx[ii].dim());
        }
        PL(nelement); PL(k); PL(i[k]); PR(_idx[k].dim()); } 
      assert(i[k] < _idx[k].dim());
    }
    int ind = i[0];
    for (int k(1); k<i.size(); k++) {
      ind *= _idx[k].dim();
      ind += i[k];
    }
    return ind;
  }

};


#endif // MULTIARRAY_HPP
