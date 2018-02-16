// Copyright (C) 2018, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license.
// See the LICENSE.txt file for details.

#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

//  Coordinate system
//
//       y
//       |
//       |
//       |
//       0-----x
//      /
//     /
//    z
//

// Cube Corners
// Corners are voxels. Numbers correspond to Morton codes of corner coordinates.
// Each cube is associated with an 8 bit mask. Each corner is assigned the bit
// at the position of its Morton code value.
//
//       2-------------------3
//      /|                  /|
//     / |                 / |
//    /  |                /  |
//   6-------------------7   |
//   |   |               |   |
//   |   |               |   |
//   |   |               |   |
//   |   |               |   |
//   |   0---------------|---1
//   |  /                |  /
//   | /                 | /
//   |/                  |/
//   4-------------------5
//


//         Cube Edges
//  
//       o--------4----------o
//      /|                  /|
//     7 |                 5 |
//    /  |                /  |
//   o--------6----------o   |
//   |   8               |   9
//   |   |               |   |
//   |   |               |   |
//   11  |               10  |
//   |   o--------0------|---o
//   |  /                |  /
//   | 3                 | 1
//   |/                  |/
//   o--------2----------o
//

/// Class which represents a cube corner of the unit cube [0,1]^3 as the
/// Morton code of its coordinates. This allows to easily get the Morton codes
/// of neighboring corners with simple XOR operations.
class CubeCornerCode {
public:
    /// Constructor initializes the code. Should be a value in {0,...,7}.
    CubeCornerCode(uint32_t c) : code(c) {}
    
    /// Get the raw Morton code of this corner.
    uint32_t getCode() const { return code;}
    
    /// Get the mask of this corner. Corresponds to the bit at the position of
    /// its Morton code value.
    uint32_t getMask() const { return 1u << code;}
    
    /// Get the neighboring code in x direction.
    uint32_t nX() const { return code ^ NEIGHBOR_X_MASK;}
    
    /// Get the neighboring code in y direction.
    uint32_t nY() const { return code ^ NEIGHBOR_Y_MASK;}
    
    /// Get the neighboring code in z direction.
    uint32_t nZ() const { return code ^ NEIGHBOR_Z_MASK;}
    
private:
    /// The Morton code
    uint32_t code;
    
    /// XOR mask for retrieving the neighbor in x direction
    static uint32_t constexpr NEIGHBOR_X_MASK = 1u;
    /// XOR mask for retrieving the neighbor in y direction
    static uint32_t constexpr NEIGHBOR_Y_MASK = 2u;
    /// XOR mask for retrieving the neighbor in z direction
    static uint32_t constexpr NEIGHBOR_Z_MASK = 4u;
};

// edge codes from the DualMC class in dualmc.h
enum class DMCEdgeCode : uint32_t {
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
    EDGE11 = 1 << 11
};

// for a corner id given by its Morton code this table gives the edge masks
// of adjacent edges in x, y, and z direction.
uint32_t cornerEdges[8][3] = {
    // {x,y,z}
    {static_cast<uint32_t>(DMCEdgeCode::EDGE0),static_cast<uint32_t>(DMCEdgeCode::EDGE8),static_cast<uint32_t>(DMCEdgeCode::EDGE3)},  // corner 0
    {static_cast<uint32_t>(DMCEdgeCode::EDGE0),static_cast<uint32_t>(DMCEdgeCode::EDGE9),static_cast<uint32_t>(DMCEdgeCode::EDGE1)},  // corner 1
    {static_cast<uint32_t>(DMCEdgeCode::EDGE4),static_cast<uint32_t>(DMCEdgeCode::EDGE8),static_cast<uint32_t>(DMCEdgeCode::EDGE7)},  // corner 2
    {static_cast<uint32_t>(DMCEdgeCode::EDGE4),static_cast<uint32_t>(DMCEdgeCode::EDGE9),static_cast<uint32_t>(DMCEdgeCode::EDGE5)},  // corner 3
    {static_cast<uint32_t>(DMCEdgeCode::EDGE2),static_cast<uint32_t>(DMCEdgeCode::EDGE11),static_cast<uint32_t>(DMCEdgeCode::EDGE3)}, // corner 4
    {static_cast<uint32_t>(DMCEdgeCode::EDGE2),static_cast<uint32_t>(DMCEdgeCode::EDGE10),static_cast<uint32_t>(DMCEdgeCode::EDGE1)}, // corner 5
    {static_cast<uint32_t>(DMCEdgeCode::EDGE6),static_cast<uint32_t>(DMCEdgeCode::EDGE11),static_cast<uint32_t>(DMCEdgeCode::EDGE7)}, // corner 6
    {static_cast<uint32_t>(DMCEdgeCode::EDGE6),static_cast<uint32_t>(DMCEdgeCode::EDGE10),static_cast<uint32_t>(DMCEdgeCode::EDGE5)}  // corner 7
};

// vector which will contain the table
std::vector<uint32_t> dualPointsList;

// Function for generating the dual marching cubes table. For each cube
// configuration it uses each corner that is clssified as inside as the starting
// corner for finding connected inside corners, that can be reached by
// traversing the cube edges. For each corner in such a connected subgraph
// we collect all edges, which connect to an outside corner.
// There is one class of configurations (126,189,219, and 231), for which this
// approach merges two original marching cubes patches into one patch.
// Example instance of this problematic class (inside 1, outside 0):
//    1------------0
//   /|           /|
//  1------------1 |
//  | |          | |
//  | |          | |
//  | |          | |
//  | 1----------|-1
//  |/           |/
//  0------------1
// Luckily, the correct patches are identical to the results of the inverted
// configurations, which are handled correctly.
void generateDualMCTable() {
    std::cout << "Generating DualMC table" << std::endl;

    // allocate space for the up to four dual point edge masks for all 256
    // configurations
    dualPointsList.resize(256 * 4);
    // corner stack for finding connected corners
    std::vector<CubeCornerCode> cornerStack;
    cornerStack.reserve(8);
    // iterate all 256 in/out cube corner cases
    for(uint32_t i = 0; i < 256; ++i) {
        // get table row of the current configuration
        uint32_t * cubeDualPointCodes = &dualPointsList[i * 4];
        
        // cube masks 0 and 255 have no intersection edges.
        // Zero row entries and continue with the next iteration.
        if(i == 0 || i == 255) {
            dualPointsList[0] = dualPointsList[1] = dualPointsList[2] = dualPointsList[3] = 0;
            continue;
        }
        
        // get the current cube configuration mask and replace the problematic
        // configurations by their inverse
        // masks
        uint32_t cubeMask = i;
        if(i == 126 || i == 189 || i == 219 || i == 231) {
            cubeMask ^= 0x000000ffu;
        }
        
        // keep track of already visited corners with a corners mask
        uint32_t processedCornersMask = 0;
        uint32_t numDualPoints = 0;
        
        // iterate all corners as start corners for finding connected corners
        for(uint32_t c = 0; c < 8; ++c) {
            CubeCornerCode startCorner(c);
            // if the corner has already been visited by a previous iteration or
            // is classified as outside we skip this iteration
            if( (startCorner.getMask() & processedCornersMask) != 0 || (cubeMask & startCorner.getMask()) == 0) {
                // mark this corner as visited
                processedCornersMask |= startCorner.getMask();
                continue;
            }
            
            // find connected corners and determin edges with surface
            // intersections. For this we initialize the subgraph traversal
            // stack with the starting corner.            
            cornerStack.emplace_back(startCorner);
            // mask for keeping track of already visited corners of the current
            // subgraph. Initialized with the starting corner
            uint32_t connectedCornersMask = startCorner.getMask();
            // variable for collecting all intersection edges of the connected
            // subgraph
            uint32_t dualPointCode = 0;
            
            // expand the connected subgraph as long as there are corners on the
            // traversal stack.
            while(!cornerStack.empty()) {
                // pop connected corner from stack
                CubeCornerCode corner = cornerStack.back();
                cornerStack.pop_back();
                assert((connectedCornersMask & corner.getMask()) != 0);
                // examine the three neighboring corners of the current corner
                std::array<CubeCornerCode,3> neighbors = {corner.nX(),corner.nY(),corner.nZ()};
                for(size_t n = 0; n < neighbors.size(); ++n) {
                    CubeCornerCode const & neighbor = neighbors[n];
                    // check if this neighboring corner is inside or outside
                    if((cubeMask & neighbor.getMask()) == 0) {
                        // have an intersection edge -> register it
                        dualPointCode |= cornerEdges[corner.getCode()][n];
                    } else {
                        // Have a connected corner. Add it to the corners stack
                        // if not already done.
                        if((connectedCornersMask & neighbor.getMask()) == 0) {
                            // register the new corner in the connected corners
                            // mask and push it on the traversal stack.
                            connectedCornersMask |= neighbor.getMask();
                            cornerStack.emplace_back(neighbor);
                        }
                    }
                }
            }
            assert((processedCornersMask & connectedCornersMask) == 0);
            // add all subgraph corners to the mask of all processed corners
            processedCornersMask |= connectedCornersMask;
            
            assert(dualPointCode != 0);
            
            // store the point code
            cubeDualPointCodes[numDualPoints] = dualPointCode;
            ++numDualPoints;
        }
        
        // set remaining point codes to zero
        for(uint32_t d = numDualPoints; d < 4; ++d) {
            cubeDualPointCodes[d] = 0;
        }
    } // END: cube code loop
}

//------------------------------------------------------------------------------

// function writing the dual marching cubes table file
void writeDualMCTable() {    
    // now write the table to a file
    char const * const filename = "generatedtable.tpp";
    std::ofstream file(filename);
    
    std::cout << "Writting DualMC table to '" << filename << '\'' << std::endl;
    
    file << "template<class T>\nint32_t DualMC<T>::dualPointsList[256][4] = {\n";
    // iterate the cube cases
    for(uint32_t cube = 0; cube < 256; ++cube) {
        uint32_t const * const codes = &dualPointsList[cube * 4];
        file << '{';
        for(int i = 0; i < 4; ++i) {
            if(i > 0)
                file << ", ";
            uint32_t code = codes[i];
            if(code != 0) {
                // extract edge codes
                bool wroteFirstEdge = false;
                int edge = 0;
                for(; code != 0; code >>= 1, ++edge) {
                    if((code & 1) == 1) {
                        if(wroteFirstEdge)
                            file << '|';
                        file << "EDGE" << edge;
                        wroteFirstEdge = true;
                    }
                }
            } else {
                file << '0';
            }
        }
        if(cube < 255)
            file << "},";
        else
            file << "}";
        file << " // " << cube << '\n';
    }
    file << "};\n";
    file.close();
}

//------------------------------------------------------------------------------

int main(int argc, char const * argv[]) {
    generateDualMCTable();
    writeDualMCTable();
    return 0;
}