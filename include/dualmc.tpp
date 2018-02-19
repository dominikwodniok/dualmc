// Copyright (C) 2017, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license. See the LICENSE.txt file for details.

/// \file   dualmc.tpp
/// \author Dominik Wodniok
/// \date   2009

//------------------------------------------------------------------------------

template<class T> inline
int DualMC<T>::getCellCode(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso) const {
    // determine for each cube corner if it is outside or inside
    int code = 0;
    if(data[gA(cx,cy,cz)] >= iso)
        code |= 1;
    if(data[gA(cx+1,cy,cz)] >= iso)
        code |= 2;
    if(data[gA(cx,cy+1,cz)] >= iso)
        code |= 4;
    if(data[gA(cx+1,cy+1,cz)] >= iso)
        code |= 8;
    if(data[gA(cx,cy,cz+1)] >= iso)
        code |= 16;
    if(data[gA(cx+1,cy,cz+1)] >= iso)
        code |= 32;
    if(data[gA(cx,cy+1,cz+1)] >= iso)
        code |= 64;
    if(data[gA(cx+1,cy+1,cz+1)] >= iso)
        code |= 128;
    return code;
}

//------------------------------------------------------------------------------

template<class T> inline
int DualMC<T>::getDualPointCode(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso, DMCEdgeCode const edge) const {
    int cubeCode = getCellCode(cx, cy, cz, iso);
    
    // is manifold dual marching cubes desired?
    if(generateManifold) {
        // The Manifold Dual Marching Cubes approach from Rephael Wenger as described in
        // chapter 3.3.5 of his book "Isosurfaces: Geometry, Topology, and Algorithms"
        // is implemente here.
        // If a problematic C16 or C19 configuration shares the ambiguous face 
        // with another C16 or C19 configuration we simply invert the cube code
        // before looking up dual points. Doing this for these pairs ensures
        // manifold meshes.
        // But this removes the dualism to marching cubes.
        
        // check if we have a potentially problematic configuration
        uint8_t const direction = problematicConfigs[uint8_t(cubeCode)];
        // If the direction code is in {0,...,5} we have a C16 or C19 configuration.
        if(direction != 255) {
            // We have to check the neighboring cube, which shares the ambiguous
            // face. For this we decode the direction. This could also be done
            // with another lookup table.
            // copy current cube coordinates into an array.
            int32_t neighborCoords[] = {cx,cy,cz};
            // get the dimension of the non-zero coordinate axis
            unsigned int const component = direction >> 1;
            // get the sign of the direction
            int32_t delta = (direction & 1) == 1 ? 1 : -1;
            // modify the correspong cube coordinate
            neighborCoords[component] += delta;
            // have we left the volume in this direction?
            if(neighborCoords[component] >= 0 && neighborCoords[component] < (dims[component]-1)) {
                // get the cube configuration of the relevant neighbor
                int neighborCubeCode = getCellCode(neighborCoords[0], neighborCoords[1], neighborCoords[2], iso);
                // Look up the neighbor configuration ambiguous face direction.
                // If the direction is valid we have a C16 or C19 neighbor.
                // As C16 and C19 have exactly one ambiguous face this face is
                // guaranteed to be shared for the pair.
                if(problematicConfigs[uint8_t(neighborCubeCode)] != 255) {
                    // replace the cube configuration with its inverse.
                    cubeCode ^= 0xff;
                }
            }
        }
    }
    for(int i = 0; i < 4; ++i)
        if(dualPointsList[cubeCode][i] & edge) {
            return dualPointsList[cubeCode][i];
        }
    return 0;
}

//------------------------------------------------------------------------------

template<class T> inline
void DualMC<T>::calculateDualPoint(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso, int const pointCode, Vertex & v) const {
    // initialize the point with lower voxel coordinates
    v.x = cx;
    v.y = cy;
    v.z = cz;
    
    // compute the dual point as the mean of the face vertices belonging to the
    // original marching cubes face
    Vertex p;
    p.x=0;
    p.y=0;
    p.z=0;
    int points = 0;

    // sum edge intersection vertices using the point code
    if(pointCode & EDGE0) {
        p.x += ((float)iso - (float)data[gA(cx,cy,cz)])/((float)data[gA(cx+1,cy,cz)]-(float)data[gA(cx,cy,cz)]);
        points++;
    }

    if(pointCode & EDGE1) {
        p.x += 1.0f;
        p.z += ((float)iso - (float)data[gA(cx+1,cy,cz)])/((float)data[gA(cx+1,cy,cz+1)]-(float)data[gA(cx+1,cy,cz)]);
        points++;
    }

    if(pointCode & EDGE2) {
        p.x += ((float)iso - (float)data[gA(cx,cy,cz+1)])/((float)data[gA(cx+1,cy,cz+1)]-(float)data[gA(cx,cy,cz+1)]);
        p.z += 1.0f;
        points++;
    }

    if(pointCode & EDGE3) {
        p.z += ((float)iso - (float)data[gA(cx,cy,cz)])/((float)data[gA(cx,cy,cz+1)]-(float)data[gA(cx,cy,cz)]);
        points++;
    }

    if(pointCode & EDGE4) {
        p.x += ((float)iso - (float)data[gA(cx,cy+1,cz)])/((float)data[gA(cx+1,cy+1,cz)]-(float)data[gA(cx,cy+1,cz)]);
        p.y += 1.0f;
        points++;
    }

    if(pointCode & EDGE5) {
        p.x += 1.0f;
        p.z += ((float)iso - (float)data[gA(cx+1,cy+1,cz)])/((float)data[gA(cx+1,cy+1,cz+1)]-(float)data[gA(cx+1,cy+1,cz)]);
        p.y += 1.0f;
        points++;
    }

    if(pointCode & EDGE6) {
        p.x += ((float)iso - (float)data[gA(cx,cy+1,cz+1)])/((float)data[gA(cx+1,cy+1,cz+1)]-(float)data[gA(cx,cy+1,cz+1)]);
        p.z += 1.0f;
        p.y += 1.0f;
        points++;
    }

    if(pointCode & EDGE7) {
        p.z += ((float)iso - (float)data[gA(cx,cy+1,cz)])/((float)data[gA(cx,cy+1,cz+1)]-(float)data[gA(cx,cy+1,cz)]);
        p.y += 1.0f;
        points++;
    }

    if(pointCode & EDGE8) {
        p.y += ((float)iso - (float)data[gA(cx,cy,cz)])/((float)data[gA(cx,cy+1,cz)]-(float)data[gA(cx,cy,cz)]);
        points++;
    }

    if(pointCode & EDGE9) {
        p.x += 1.0f;
        p.y += ((float)iso - (float)data[gA(cx+1,cy,cz)])/((float)data[gA(cx+1,cy+1,cz)]-(float)data[gA(cx+1,cy,cz)]);
        points++;
    }

    if(pointCode & EDGE10) {
        p.x += 1.0f;
        p.y += ((float)iso - (float)data[gA(cx+1,cy,cz+1)])/((float)data[gA(cx+1,cy+1,cz+1)]-(float)data[gA(cx+1,cy,cz+1)]);
        p.z += 1.0f;
        points++;
    }

    if(pointCode & EDGE11) {
        p.z += 1.0f;
        p.y += ((float)iso - (float)data[gA(cx,cy,cz+1)])/((float)data[gA(cx,cy+1,cz+1)]-(float)data[gA(cx,cy,cz+1)]);
        points++;
    }

    // divide by number of accumulated points
    float invPoints = 1.0f / (float)points;
    p.x*= invPoints;
    p.y*= invPoints;
    p.z*= invPoints;

    // offset point by voxel coordinates
    v.x += p.x;
    v.y += p.y;
    v.z += p.z;
}

//------------------------------------------------------------------------------

template<class T> inline
QuadIndexType DualMC<T>::getSharedDualPointIndex(
    int32_t const cx, int32_t const cy, int32_t const cz,
    VolumeDataType const iso, DMCEdgeCode const edge,
    std::vector<Vertex> & vertices
    ) {
    // create a key for the dual point from its linearized cell ID and point code
    DualPointKey key;
    key.linearizedCellID = gA(cx,cy,cz);
    key.pointCode = getDualPointCode(cx,cy,cz,iso,edge);
    
    // have we already computed the dual point?
    auto iterator = pointToIndex.find(key);
    if(iterator != pointToIndex.end()) {
        // just return the dual point index
        return iterator->second;
    } else {
        // create new vertex and vertex id
        QuadIndexType newVertexId = vertices.size();
        vertices.emplace_back();
        calculateDualPoint(cx,cy,cz,iso,key.pointCode, vertices.back());
        // insert vertex ID into map and also return it
        pointToIndex[key] = newVertexId;
        return newVertexId;
    }
}

//------------------------------------------------------------------------------

template<class T> inline
void DualMC<T>::build(
    VolumeDataType const * data,
    int32_t const dimX, int32_t const dimY, int32_t const dimZ,
    VolumeDataType const iso,
    bool const generateManifold,
    bool const generateSoup,
    std::vector<Vertex> & vertices,
    std::vector<Quad> & quads
    ) {

    // set members
    this->dims[0] = dimX;
    this->dims[1] = dimY;
    this->dims[2] = dimZ;
    this->data = data;
    this->generateManifold = generateManifold;
    
    // clear vertices and quad indices
    vertices.clear();
    quads.clear();
    
    // generate quad soup or shared vertices quad list
    if(generateSoup) {
        buildQuadSoup(iso,vertices,quads);
    } else {
        buildSharedVerticesQuads(iso,vertices,quads);
    }
}

//------------------------------------------------------------------------------

template<class T> inline
void DualMC<T>::buildQuadSoup(
    VolumeDataType const iso,
    std::vector<Vertex> & vertices,
    std::vector<Quad> & quads
    ) {
    
    int32_t const reducedX = dims[0] - 2;
    int32_t const reducedY = dims[1] - 2;
    int32_t const reducedZ = dims[2] - 2;

    Vertex vertex0;
    Vertex vertex1;
    Vertex vertex2;
    Vertex vertex3;
    int pointCode;

    // iterate voxels
    for(int32_t z = 0; z < reducedZ; ++z)
        for(int32_t y = 0; y < reducedY; ++y)
            for(int32_t x = 0; x < reducedX; ++x) {
                // construct quad for x edge
                if(z > 0 && y > 0) {
                    // is edge intersected?
                    bool const entering = data[gA(x,y,z)] < iso && data[gA(x+1,y,z)] >= iso;
                    bool const exiting  = data[gA(x,y,z)] >= iso && data[gA(x+1,y,z)] < iso;
                    if(entering || exiting){
                        // generate quad
                        pointCode = getDualPointCode(x,y,z,iso,EDGE0);
                        calculateDualPoint(x,y,z,iso,pointCode, vertex0);

                        pointCode = getDualPointCode(x,y,z-1,iso,EDGE2);
                        calculateDualPoint(x,y,z-1,iso,pointCode, vertex1);

                        pointCode = getDualPointCode(x,y-1,z-1,iso,EDGE6);
                        calculateDualPoint(x,y-1,z-1,iso,pointCode, vertex2);

                        pointCode = getDualPointCode(x,y-1,z,iso,EDGE4);
                        calculateDualPoint(x,y-1,z,iso,pointCode, vertex3);
                        
                        if(entering) {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex1);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex3);
                        } else {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex3);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex1);
                        }
                    }
                }
                
                // construct quad for y edge
                if(z > 0 && x > 0) {
                    // is edge intersected?
                    bool const entering = data[gA(x,y,z)] < iso && data[gA(x,y+1,z)] >= iso;
                    bool const exiting  = data[gA(x,y,z)] >= iso && data[gA(x,y+1,z)] < iso;
                    if(entering || exiting){
                        // generate quad
                        pointCode = getDualPointCode(x,y,z,iso,EDGE8);
                        calculateDualPoint(x,y,z,iso,pointCode, vertex0);

                        pointCode = getDualPointCode(x,y,z-1,iso,EDGE11);
                        calculateDualPoint(x,y,z-1,iso,pointCode, vertex1);

                        pointCode = getDualPointCode(x-1,y,z-1,iso,EDGE10);
                        calculateDualPoint(x-1,y,z-1,iso,pointCode, vertex2);

                        pointCode = getDualPointCode(x-1,y,z,iso,EDGE9);
                        calculateDualPoint(x-1,y,z,iso,pointCode, vertex3);
                        
                        if(exiting) {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex1);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex3);
                        } else {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex3);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex1);
                        }
                    }
                }

                // construct quad for z edge
                if(x > 0 && y > 0) {
                    // is edge intersected?
                    bool const entering = data[gA(x,y,z)] < iso && data[gA(x,y,z+1)] >= iso;
                    bool const exiting  = data[gA(x,y,z)] >= iso && data[gA(x,y,z+1)] < iso;
                    if(entering || exiting){
                        // generate quad
                        pointCode = getDualPointCode(x,y,z,iso,EDGE3);
                        calculateDualPoint(x,y,z,iso,pointCode, vertex0);
                        
                        pointCode = getDualPointCode(x-1,y,z,iso,EDGE1);
                        calculateDualPoint(x-1,y,z,iso,pointCode, vertex1);

                        pointCode = getDualPointCode(x-1,y-1,z,iso,EDGE5);
                        calculateDualPoint(x-1,y-1,z,iso,pointCode, vertex2);

                        pointCode = getDualPointCode(x,y-1,z,iso,EDGE7);
                        calculateDualPoint(x,y-1,z,iso,pointCode, vertex3);
                        
                        if(exiting) {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex1);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex3);
                        } else {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex3);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex1);
                        }
                    }
                }
            }
    
    // generate triangle soup quads
    size_t const numQuads = vertices.size() / 4;
    quads.reserve(numQuads);
    for (size_t i = 0; i < numQuads; ++i) {
        quads.emplace_back(i * 4, i * 4 + 1, i * 4 + 2, i * 4 + 3);
    }
}

//------------------------------------------------------------------------------

template<class T> inline
void DualMC<T>::buildSharedVerticesQuads(
    VolumeDataType const iso,
    std::vector<Vertex> & vertices,
    std::vector<Quad> & quads
    ) {
    

    int32_t const reducedX = dims[0] - 2;
    int32_t const reducedY = dims[1] - 2;
    int32_t const reducedZ = dims[2] - 2;

    QuadIndexType i0,i1,i2,i3;
    
    pointToIndex.clear();

    // iterate voxels
    for(int32_t z = 0; z < reducedZ; ++z)
        for(int32_t y = 0; y < reducedY; ++y)
            for(int32_t x = 0; x < reducedX; ++x) {
                // construct quads for x edge
                if(z > 0 && y > 0) {
                    bool const entering = data[gA(x,y,z)] < iso && data[gA(x+1,y,z)] >= iso;
                    bool const exiting  = data[gA(x,y,z)] >= iso && data[gA(x+1,y,z)] < iso;
                    if(entering || exiting){
                        // generate quad
                        i0 = getSharedDualPointIndex(x,y,z,iso,EDGE0,vertices);
                        i1 = getSharedDualPointIndex(x,y,z-1,iso,EDGE2,vertices);
                        i2 = getSharedDualPointIndex(x,y-1,z-1,iso,EDGE6,vertices);
                        i3 = getSharedDualPointIndex(x,y-1,z,iso,EDGE4,vertices);
                        
                        if(entering) {
                            quads.emplace_back(i0,i1,i2,i3);
                        } else {
                            quads.emplace_back(i0,i3,i2,i1);
                        }
                    }
                }
                
                // construct quads for y edge
                if(z > 0 && x > 0) {
                    bool const entering = data[gA(x,y,z)] < iso && data[gA(x,y+1,z)] >= iso;
                    bool const exiting  = data[gA(x,y,z)] >= iso && data[gA(x,y+1,z)] < iso;
                    if(entering || exiting){
                        // generate quad
                        i0 = getSharedDualPointIndex(x,y,z,iso,EDGE8,vertices);
                        i1 = getSharedDualPointIndex(x,y,z-1,iso,EDGE11,vertices);
                        i2 = getSharedDualPointIndex(x-1,y,z-1,iso,EDGE10,vertices);
                        i3 = getSharedDualPointIndex(x-1,y,z,iso,EDGE9,vertices);
                        
                        if(exiting) {
                            quads.emplace_back(i0,i1,i2,i3);
                        } else {
                            quads.emplace_back(i0,i3,i2,i1);
                        }
                    }
                }

                // construct quads for z edge
                if(x > 0 && y > 0) {
                    bool const entering = data[gA(x,y,z)] < iso && data[gA(x,y,z+1)] >= iso;
                    bool const exiting  = data[gA(x,y,z)] >= iso && data[gA(x,y,z+1)] < iso;
                    if(entering || exiting){
                        // generate quad
                        i0 = getSharedDualPointIndex(x,y,z,iso,EDGE3,vertices);                        
                        i1 = getSharedDualPointIndex(x-1,y,z,iso,EDGE1,vertices);
                        i2 = getSharedDualPointIndex(x-1,y-1,z,iso,EDGE5,vertices);
                        i3 = getSharedDualPointIndex(x,y-1,z,iso,EDGE7,vertices);
                        
                        if(exiting) {
                            quads.emplace_back(i0,i1,i2,i3);
                        } else {
                            quads.emplace_back(i0,i3,i2,i1);
                        }
                    }
                } 
            }
}
