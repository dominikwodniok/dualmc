// Copyright (C) 2017, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license.
// See the LICENSE.txt file for details.

#ifndef DUALMC_H_INCLUDED
#define DUALMC_H_INCLUDED

/// \file   dualmc.h
/// \author Dominik Wodniok
/// \date   2009

// c includes
#include <cstdint>

// stl includes
#include <unordered_map>
#include <vector>

namespace dualmc {
    

typedef float VertexComponentsType;
typedef int32_t QuadIndexType;

/// vertex structure for dual points
struct Vertex {
    /// non-initializing constructor
    Vertex();

    /// initializing constructor
    Vertex(VertexComponentsType x, VertexComponentsType y, VertexComponentsType z);
    
    /// initializing constructor
    Vertex(Vertex const & v);
    
    // components
    VertexComponentsType x,y,z;
};

/// quad indices structure
struct Quad {
    /// non-initializing constructor
    Quad();
    
    /// initializing constructor
    Quad(QuadIndexType i0, QuadIndexType i1,QuadIndexType i2, QuadIndexType i3);
    
    // quad indices
    QuadIndexType i0,i1,i2,i3;
    
};

/// \class  DualMC
/// \author Dominik Wodniok
/// \date   2009
/// Class which implements the dual marching cubes algorithm from Gregory M. Nielson.
/// Faces and vertices of the standard marching cubes algorithm correspond to
/// vertices and faces in the dual algorithm. As a vertex in standard marching cubes
/// usually is shared by 4 faces, the dual mesh is entirely made from quadrangles.
/// Unfortunately, under rare circumstances the original algorithm can create
/// non-manifold meshes. See the remarks of the original paper on this.
/// The class optionally can guarantee manifold meshes by taking the Manifold
/// Dual Marching Cubes approach from Rephael Wenger as described in
/// chapter 3.3.5 of his book "Isosurfaces: Geometry, Topology, and Algorithms".
template<class T> class DualMC {
public:
    // typedefs
    typedef T VolumeDataType;

    /// Extracts the iso surface for a given volume and iso value.
    /// Output is a list of vertices and a list of indices, which connect
    /// vertices to quads.
    /// The quad mesh either uses shared vertex indices or is a quad soup if
    /// desired.
    void build(
        VolumeDataType const * data,
        int32_t const dimX, int32_t const dimY, int32_t const dimZ,
        VolumeDataType const iso,
        bool const generateManifold,
        bool const generateSoup,
        std::vector<Vertex> & vertices,
        std::vector<Quad> & quads
        );

private:

    /// Extract quad mesh with shared vertex indices.
    void buildSharedVerticesQuads(
        VolumeDataType const iso,
        std::vector<Vertex> & vertices,
        std::vector<Quad> & quads
        );
        
    /// Extract quad soup.
    void buildQuadSoup(
        VolumeDataType const iso,
        std::vector<Vertex> & vertices,
        std::vector<Quad> & quads
        );


private:

    /// enum with edge codes for a 12-bit voxel edge mask to indicate
    /// grid edges which intersect the ISO surface of classic marching cubes
    enum DMCEdgeCode {
        EDGE0 = 1,
        EDGE1 = 1 << 1,
        EDGE2 = 1 << 2,
        EDGE3 = 1 << 3,
        EDGE4 = 1 << 4,
        EDGE5 = 1 << 5,
        EDGE6 = 1 << 6,
        EDGE7 = 1 << 7,
        EDGE8 = 1 << 8,
        EDGE9 = 1 << 9,
        EDGE10 = 1 << 10,
        EDGE11 = 1 << 11,
        FORCE_32BIT = 0xffffffff
    };

    /// get the 8-bit in-out mask for the voxel corners of the cell cube at (cx,cy,cz)
    /// and the given iso value
    int getCellCode(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso) const;

    /// Get the 12-bit dual point code mask, which encodes the traditional
    /// marching cube vertices of the traditional marching cubes face which
    /// corresponds to the dual point.
    /// This is also where the manifold dual marching cubes algorithm is
    /// implemented.
    int getDualPointCode(int32_t const cx, int32_t const cy, int32_t const cz,
      VolumeDataType const iso, DMCEdgeCode const edge) const;

    /// Given a dual point code and iso value, compute the dual point.
    void calculateDualPoint(int32_t const cx, int32_t const cy, int32_t const cz,
      VolumeDataType const iso, int const pointCode, Vertex &v) const;

    /// Get the shared index of a dual point which is uniquly identified by its
    /// cell cube index and a cube edge. The dual point is computed,
    /// if it has not been computed before.
    QuadIndexType getSharedDualPointIndex(int32_t const cx, int32_t const cy, int32_t const cz,
      VolumeDataType const iso, DMCEdgeCode const edge,
      std::vector<Vertex> & vertices);
    
    /// Compute a linearized cell cube index.
    int32_t gA(int32_t const x, int32_t const y, int32_t const z) const;

private:
    // static lookup tables needed for (manifold) dual marching cubes

    /// Dual Marching Cubes table
    /// Encodes the edge vertices for the 256 marching cubes cases.
    /// A marching cube case produces up to four faces and ,thus, up to four
    /// dual points.
    static int32_t const dualPointsList[256][4];
    
    /// Table which encodes the ambiguous face of cube configurations, which
    /// can cause non-manifold meshes.
    /// Needed for manifold dual marching cubes.
    static uint8_t const problematicConfigs[256];
    
private:

    /// convenience volume extent array for x-,y-, and z-dimension
    int32_t dims[3];

    /// convenience volume data point
    VolumeDataType const * data;
    
    /// store whether the manifold dual marching cubes algorithm should be
    /// applied.
    bool generateManifold;
    
    /// Dual point key structure for hashing of shared vertices
    struct DualPointKey {
        // a dual point can be uniquely identified by ite linearized volume cell
        // id and point code
        int32_t linearizedCellID;
        int pointCode;
        /// Equal operator for unordered map
        bool operator==(DualPointKey const & other) const;
    };
    
    /// Functor for dual point key hash generation
    struct DualPointKeyHash {
        size_t operator()(DualPointKey const & k) const {
            return size_t(k.linearizedCellID) | (size_t(k.pointCode) << 32u);
        }
    };
    
    /// Hash map for shared vertex index computations
    std::unordered_map<DualPointKey,QuadIndexType,DualPointKeyHash> pointToIndex;
};

// inline function definitions

//------------------------------------------------------------------------------

inline
Vertex::Vertex(){}

//------------------------------------------------------------------------------

inline
Vertex::Vertex(
    VertexComponentsType x,
    VertexComponentsType y,
    VertexComponentsType z
    ) : x(x), y(y), z(z) {}

//------------------------------------------------------------------------------

inline
Vertex::Vertex(Vertex const & v) : x(v.x), y(v.y), z(v.z) {}

//------------------------------------------------------------------------------

inline
Quad::Quad(){}

//------------------------------------------------------------------------------

inline
Quad::Quad(
    QuadIndexType i0,
    QuadIndexType i1,
    QuadIndexType i2,
    QuadIndexType i3
    ) : i0(i0),i1(i1),i2(i2),i3(i3) {}

//------------------------------------------------------------------------------

template<class T> inline
int32_t DualMC<T>::gA(int32_t const x, int32_t const y, int32_t const z) const {
    return x + dims[0] * (y + dims[1] * z);
}

//------------------------------------------------------------------------------
template<class T> inline
bool DualMC<T>::DualPointKey::operator==(typename DualMC<T>::DualPointKey const & other) const {
    return linearizedCellID == other.linearizedCellID && pointCode == other.pointCode;
}

#include "dualmc.tpp"

#include "dualmc_tables.tpp"

} // END: namespace dualmc
#endif // DUALMC_H_INCLUDED
