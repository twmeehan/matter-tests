#ifndef CSV2ABC_HPP
#define CSV2ABC_HPP

#include <cmath>
#include <vector>

#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreOgawa/All.h>
#include <ImathRandom.h>
#include <Alembic/AbcCoreAbstract/Tests/Assert.h>

namespace AbcG = Alembic::AbcGeom;
using namespace AbcG;

// See Simulation Replay!!!!

void csv2abc(/*std::string csv_filename, */ std::string abc_filename){

        int numParticles = 5;

        std::vector<Alembic::Util::uint64_t> my_ids;
        std::vector<V3f> my_positions;
        std::vector<V3f> my_velocities;

        for (int i = 0; i < numParticles; ++i){
            my_positions.push_back(V3f(i, i, i));
            my_velocities.push_back(V3f(100.0-i, 0, 0));
            my_ids.push_back(i);
        }

        for (int i = 0; i < numParticles; ++i)
            std::cout << my_positions[i][0] << " " << my_positions[i][1] << " " << my_positions[i][2]  << std::endl;

        OArchive archive( Alembic::AbcCoreOgawa::WriteArchive(), abc_filename );
        OObject topObj( archive, kTop );
        size_t iNumFrames = 20;
        chrono_t iFps = 1.0/24.0;

        TimeSampling ts(iFps, 0.0);
        Alembic::Util::uint32_t tsidx = topObj.getArchive().addTimeSampling(ts);

        OPoints partsOut( topObj, "simpleParticles", tsidx );
        std::cout << "Created Simple Particles" << std::endl;

        OPointsSchema pSchema = partsOut.getSchema();
        MetaData mdata;
        SetGeometryScope( mdata, kVaryingScope );
        OV3fArrayProperty velOut( pSchema, "velocity", mdata, tsidx );

        // Get seconds per frame.
        chrono_t iSpf = 1.0 / iFps;

        // Loop over the frames.
        for ( index_t sampIndex = 0; sampIndex < (index_t)iNumFrames; ++sampIndex )
        {
            // First, write the sample.
            // OPointsSchema::Sample psamp(
            //     V3fArraySample( my_positions ),
            //     UInt64ArraySample( my_ids ) );
            // pSchema.set( psamp );
            velOut.set( V3fArraySample( my_velocities ) );

            // // Now time step.
            // parts.timeStep( iSpf );

            std::cout << "Wrote to frame: " << sampIndex << std::endl;
        }


} // end csv2abc



#endif  // CSV2ABC_HPP
