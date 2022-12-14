#ifndef INC_FRAMEARRAY_H
#define INC_FRAMEARRAY_H
#include "Frame.h"
/// Hold a basic array of Frames
class FrameArray {
  public:
    FrameArray() {}
    void resize(int nIn)                   { farray_.resize(nIn);     }
    Frame&       operator[](int idx)       { return farray_[idx];     }
    Frame const& operator[](int idx) const { return farray_[idx];     }
    void AddFrame( const Frame& fIn )      { farray_.push_back( fIn );}
    size_t Size()                    const { return farray_.size();   }

    int SetupFrames(std::vector<Atom> const& Atoms, CoordinateInfo const& cInfo) {
      for (std::vector<Frame>::iterator myF = farray_.begin(); myF != farray_.end(); ++myF)
        if (myF->SetupFrameV(Atoms, cInfo) != 0) return 1;
      return 0;
    }
 
    typedef std::vector<Frame>::iterator iterator;
    iterator begin() { return farray_.begin(); }
    iterator end()   { return farray_.end();   }
  private:
    std::vector<Frame> farray_;
};
#endif
