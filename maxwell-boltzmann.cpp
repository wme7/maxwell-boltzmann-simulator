/* +---------------------------------------------------------------------------+
   |          Maxwell-Boltzmann distribution demonstration (2D/3D)             |
   |                                                                           |
   |   Copyright (C) 2012  Jose Luis Blanco                                    |
   |                                                                           |
   |    This software was written by Jose Luis Blanco, University of Malaga    |
   |       for publication in http://www.ciencia-explicada.com/                |
   |    Contact: blog.ciencia.explicada@gmail.com                              |
   |                                                                           |
   |  This program is free software: you can redistribute it and/or modify     |
   |     it under the terms of the GNU General Public License as published by  |
   |     the Free Software Foundation, either version 3 of the License, or     |
   |     (at your option) any later version.                                   |
   |                                                                           |
   |  This program is distributed in the hope that it will be useful,          |
   |     but WITHOUT ANY WARRANTY; without even the implied warranty of        |
   |     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
   |     GNU General Public License for more details.                          |
   |                                                                           |
   |     You should have received a copy of the GNU General Public License     |
   |     along with this.  If not, see <http://www.gnu.org/licenses/>.         |
   |                                                                           |
   +---------------------------------------------------------------------------+ */

#include <mrpt/base.h>
#include <mrpt/opengl.h>
#include <mrpt/gui.h>

using namespace mrpt;
using namespace std;
using namespace mrpt::gui;
using namespace mrpt::random;
using namespace mrpt::opengl;
using namespace mrpt::math;
using namespace mrpt::utils;

// This is the key macro to switch 2D <-> 3D
// Note: it will be defined via a compiler flag from CMake, don't uncomment.
//#define  SIMUL_2D


#ifdef SIMUL_2D
	const size_t  	N_MASSES = 50;
	const size_t	DOFs     = 2;
#else
	const size_t  	N_MASSES = 1500;
	const size_t	DOFs     = 3;
#endif

const double	BOX      = 1;  // Simulation box size (meters)
const double 	RADIUS   = 10e-3;  // Particle radius (meters)

const double   MASS_amu = 1.6605389e-27; // kg
const double   MASS    = 16*MASS_amu;  // Approx: Nitrogen. Only to evaluate the temperature

const size_t nBins = 50;  // To draw the histogram
const int DEFAULT_WIN_WIDTH = 600;

const double INIT_VEL = 1.0;

struct TMass
{
	TMass() : x(0),y(0),z(0),vx(0),vy(0),vz(0)
	{ }

	double	x,y,z;
	double	vx,vy,vz;
	double  radius;
#ifdef 	SIMUL_2D
	opengl::CDiskPtr	obj3d;
#endif

};


typedef vector<TMass> TMassVector;

template <typename Derived>
struct ParticleCloudAdaptor
{
	typedef double coord_t;

	const Derived &obj; //!< A const ref to the data set origin

	/// The constructor that sets the data set source
	ParticleCloudAdaptor(const Derived &obj_) : obj(obj_) { }

	/// CRTP helper method
	inline const Derived& derived() const { return obj; }

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return derived().size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline coord_t kdtree_distance(const coord_t *p1, const size_t idx_p2,size_t size) const
	{
		const coord_t d0=p1[0]-derived()[idx_p2].x;
		const coord_t d1=p1[1]-derived()[idx_p2].y;
		const coord_t d2=p1[2]-derived()[idx_p2].z;
		return d0*d0+d1*d1+d2*d2;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline coord_t kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim==0) return derived()[idx].x;
		else if (dim==1) return derived()[idx].y;
		else return derived()[idx].z;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const { return false; }
};


struct TStats
{
	TStats() :
		overall_Ec(0),
		total_presure_energy(0),
		max_vel(0),
		maxwell_vel_std(0),
		temperature(0)
	{
	}

	double 		overall_Ec; // Kinetic energy.
	double 		total_presure_energy; // Kinetic energy.
	vector_double	histogram, hist_centers;      // Histogram of velocities
	double         max_vel;
	double 		maxwell_vel_std;
	double         temperature;
};

void simulate( TMassVector &objs, double At, TStats &s);
void compute_stats( const TMassVector &objs, TStats &s);
void draw_stats( opengl::COpenGLViewportPtr & gl_extra, const TStats &stats );

// ------------------------------------------------------
//				BoltzmannDemo
// ------------------------------------------------------
void BoltzmannDemo()
{
	double SPEED_UP_FACTOR = 1;
	double LARGEST_STEP = 10e-3;
	bool    is_tracking = false; // Draw the random walk path of the first particle.

	// Usage:
	std::cout << 
	"Recognized key commands:\n"
	" 0   : Add new static particle.\n"
	" 1-4 : Launch particles of increasing speeds.\n"
	" 6   : Launch a very fast (x100) particle.\n"
	" 8-9 : Launch a lot of new particles.\n"
	" t   : Toogle tracking of random-walk.\n"
	" x/z : Simulation speed up/down (limited by CPU time).\n"
	" +/- : Increase/decrease point particle size.\n"
	"\n";



	CDisplayWindow3D	win("Boltzmann simulator",DEFAULT_WIN_WIDTH,700);

	randomGenerator.randomize();

#ifdef SIMUL_2D
	win.setCameraElevationDeg( 90.0f );
	win.setCameraAzimuthDeg( 0.0f );
	win.setCameraZoom( 3 );
	win.setCameraPointingToPoint(.5*BOX,.5*BOX,0);
#else
	win.setCameraElevationDeg( 50.0f );
	win.setCameraAzimuthDeg( 30.0f );
	win.setCameraZoom( 5 );
#endif

	COpenGLScenePtr &theScene = win.get3DSceneAndLock();

	COpenGLViewportPtr gl_view_main = theScene->getViewport("main");
	gl_view_main->setCustomBackgroundColor(TColorf(0.1f,0.1f,0.1f));
	gl_view_main->setViewportPosition(0,150,1.0,0.8);

	COpenGLViewportPtr gl_extra = theScene->createViewport("extra");
	gl_extra->setViewportPosition(0,0,1.0,150);

#ifndef 	SIMUL_2D
	opengl::CPointCloudColouredPtr gl_particles = opengl::CPointCloudColoured::Create();
	gl_particles->setPointSize(2.0);
	gl_particles->enablePointSmooth();
#endif

	// Draw (optionally) the random walk of one particle:
	mrpt::opengl::CSetOfLinesPtr gl_random_walk = mrpt::opengl::CSetOfLines::Create();
	gl_random_walk->setLineWidth(1);
	gl_random_walk->setColor_u8( TColor(220,220,220, 180) );  // RGBA

	// Modify the scene:
	// ------------------------------------------------------
	{
		//opengl::CGridPlaneXYPtr obj = opengl::CGridPlaneXY::Create(-BOX,BOX+1e-03, -BOX,BOX+1e-03,  0,0.1);
		opengl::CBoxPtr obj = opengl::CBox::Create(
			TPoint3D(0,0,0),
#ifdef 	SIMUL_2D
			TPoint3D(BOX,BOX,0)
#else
			TPoint3D(BOX,BOX,BOX)
#endif
			);
		obj->setColor(0.3,0.3,0.3);
		obj->setWireframe(true);
		obj->setLineWidth(1.0);
		theScene->insert( obj );

		theScene->insert(gl_random_walk);

#ifndef 	SIMUL_2D
		theScene->insert(gl_particles);
#endif
	}

	// Create the masses:
	// ----------------------------------------------------
	TMassVector	masses(N_MASSES);

	// Init at random poses & create opengl objects:
	for (size_t i=0;i<N_MASSES;i++)
	{
		masses[i].x = randomGenerator.drawUniform(0,BOX);
		masses[i].y = randomGenerator.drawUniform(0,BOX);
#ifdef 	SIMUL_2D
		masses[i].z = 0;
#else
		masses[i].z = randomGenerator.drawUniform(0,BOX);
#endif

		masses[i].vx = 0;
		masses[i].vy = 0;
		masses[i].vz = 0;

		masses[i].radius = RADIUS;

#ifdef 	SIMUL_2D
		opengl::CDiskPtr & obj = masses[i].obj3d = opengl::CDisk::Create(RADIUS,0, 10,1);

		obj->setColor(
			randomGenerator.drawUniform(0.1,0.9),
			randomGenerator.drawUniform(0.1,0.9),
			randomGenerator.drawUniform(0.1,0.9)  );

		obj->setLocation( masses[i].x, masses[i].y, masses[i].z );
		theScene->insert( obj );
#else
		gl_particles->push_back(
			masses[i].x, masses[i].y, masses[i].z,
			randomGenerator.drawUniform(0.3,1.0),randomGenerator.drawUniform(0.3,1.0),randomGenerator.drawUniform(0.3,1.0)
			);
#endif
	}

	// IMPORTANT!!! IF NOT UNLOCKED, THE WINDOW WILL NOT BE UPDATED!
	win.unlockAccess3DScene();


	double smooth_pressure = 0;
	const double PRESS_ALPHA = 0.999;


	mrpt::utils::CTicTac	tictac;
	tictac.Tic();

	double t0 = tictac.Tac();

	double simul_time = 0;

	while (win.isOpen() )
	{
		double t1 = tictac.Tac();
		double At = (t1-t0) * SPEED_UP_FACTOR;
		simul_time+=At;
		t0 = t1;

		// Simulate a At, possibly in several small steps:
		// ------------------------------------------------
		TStats stats;

		size_t  n_steps = ceil(At/LARGEST_STEP)+1;
		double At_steps = At / n_steps;
		n_steps = min(n_steps,size_t(3));
		for (size_t j=0;j<n_steps;j++)
			simulate( masses, At_steps,stats);


		// Update stats:
		compute_stats(masses, stats);

		// Pressure:
		const double current_pressure = stats.total_presure_energy;
		smooth_pressure = PRESS_ALPHA*smooth_pressure + (1-PRESS_ALPHA)*current_pressure;
		const size_t nParts = masses.size();

		// Lock 3D view:
		win.get3DSceneAndLock();

		// Refresh energy stats graphs:
		static int cnt_refresh = 0;
		if (cnt_refresh++>10)
		{
			cnt_refresh = 0;
			draw_stats(gl_extra, stats);
		}

		// Show random walk?
		if (is_tracking)
		{
			mrpt::math::TSegment3D  sg;
			sg.point2.x = masses[0].x;
			sg.point2.y = masses[0].y;
			sg.point2.z = masses[0].z;

			if (gl_random_walk->getLineCount()!=0)
					sg.point1 = gl_random_walk->rbegin()->point2;
			else	sg.point1 = sg.point2;

			gl_random_walk->appendLine(sg);
		}

		// Update texts:
		win.addTextMessage(
			5,-20, format("Kin.E.=%.02eJ Pressure=%.02ePa Max.Vel.=%.03fm/s", stats.overall_Ec, 1e9*smooth_pressure, stats.max_vel),
			TColorf(1,1,1),
			"mono",10, mrpt::opengl::NICE,
			100 );

		win.addTextMessage(
			5,-35, format("%u parts. Simul. time:%9.03fs  Speed: x%.4f", static_cast<unsigned int>(nParts), simul_time, SPEED_UP_FACTOR),
			TColorf(1,1,1),
			"mono",10, mrpt::opengl::NICE,
			101 );


		// Adaptative simulation timesteps:
		if (stats.max_vel>0)
			LARGEST_STEP = RADIUS * 0.3 / stats.max_vel;

		// Update particles position in the graphs:
		for (size_t i=0;i<masses.size();i++)
		{
#ifdef 	SIMUL_2D
			opengl::CDiskPtr & obj = masses[i].obj3d;
			obj->setLocation( masses[i].x, masses[i].y, masses[i].z );
#else
			gl_particles->setPoint_fast(i,masses[i].x, masses[i].y, masses[i].z);
#endif
		}

		// IMPORTANT!!! IF NOT UNLOCKED, THE WINDOW WILL NOT BE UPDATED!
		win.unlockAccess3DScene();

		// Update window:
		win.forceRepaint();
		mrpt::system::sleep(1); 

		if (win.keyHit())
		{
			const int ch = win.getPushedKey();

			switch (ch)
			{
#ifndef 	SIMUL_2D
			case '+':
				gl_particles->setPointSize( gl_particles->getPointSize() + 1 );
				break;
			case '-':
				if (gl_particles->getPointSize()>1)
					gl_particles->setPointSize( gl_particles->getPointSize() - 1 );
				break;
#endif

			case 't':
				if (is_tracking) {
					is_tracking = false;
				}
				else {
					is_tracking = true;
					gl_random_walk->clear();
				}
				break;

			case 'z':
				SPEED_UP_FACTOR*=0.5;
				break;
			case 'x':
				SPEED_UP_FACTOR*=2;
				break;

			case '0':
			case '1':
			case '2':
			case '3':
			case '4':

			case '6': // x100 speed

			case '8':
			case '9':
				{
					COpenGLScenePtr &scene = win.get3DSceneAndLock();

					const size_t nNewParts = (ch=='8' || ch=='9') ? 25 : 1;
					const double V = (ch=='8' || ch=='9') ?   INIT_VEL * 2*(ch-'7') : ( ch=='6' ? INIT_VEL * 100 : INIT_VEL * (ch-'0') );

					for (size_t k=0;k<nNewParts;k++)
					{
						masses.resize(masses.size()+1);
						const size_t i = masses.size()-1;

						masses[i].x = randomGenerator.drawUniform(0,BOX);
						masses[i].y = ch!='0' ? -BOX*0.5  : randomGenerator.drawUniform(0,BOX);
	#ifdef 	SIMUL_2D
						masses[i].z = 0;
	#else
						masses[i].z = randomGenerator.drawUniform(0,BOX);
	#endif


						masses[i].vx = V * randomGenerator.drawUniform(0,0.1);
						masses[i].vy = V;
	#ifdef 	SIMUL_2D
						masses[i].vz = 0;
	#else
						masses[i].vz = V * randomGenerator.drawUniform(0,0.1);
	#endif

						masses[i].radius = RADIUS;

	#ifdef 	SIMUL_2D
						opengl::CDiskPtr & obj = masses[i].obj3d = opengl::CDisk::Create(RADIUS,0, 10,1);
						obj->setColor(
							randomGenerator.drawUniform(0.2,0.9),
							randomGenerator.drawUniform(0.2,0.9),
							randomGenerator.drawUniform(0.2,0.9)  );

						obj->setLocation( masses[i].x, masses[i].y, masses[i].z );
						scene->insert( obj );
	#else
						gl_particles->push_back(
							masses[i].x, masses[i].y, masses[i].z,
							randomGenerator.drawUniform(0.3,1.0),randomGenerator.drawUniform(0.3,1.0),randomGenerator.drawUniform(0.3,1.0)
							);
	#endif
					}

					win.unlockAccess3DScene();
				}
				break;
			};
		}

	};
}

void simulate( TMassVector &objs, double At, TStats &s)
{
	const size_t N=objs.size();

	// Build KD-tree
	// ------------------------------------
	typedef ParticleCloudAdaptor<TMassVector> PC2KD;
	const PC2KD  pc2kd(objs); // The adaptor

	// construct a kd-tree index:
	typedef nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<double, PC2KD > ,
		PC2KD,
		3 /* dim */
		> my_kd_tree_t;

	my_kd_tree_t   index(3 /*dim*/, pc2kd, nanoflann::KDTreeSingleIndexAdaptorParams() );
	index.buildIndex();

	const TPoint3D dir_up(0,0,1);

	// Test for ellastic collisions:
	// ------------------------------------

	vector<bool> already_collide;
	already_collide.assign(N, false);

	for (size_t i=0;i<N;i++)
	{
		TMass &m = objs[i];

		// Borders?
		double press_v = 0;
		if (m.x+m.radius > BOX && m.vx>0) { m.vx=-m.vx; press_v += m.vx; }
		if (m.y+m.radius > BOX && m.vy>0) { m.vy=-m.vy; press_v += m.vy; }
		if (m.z+m.radius > BOX && m.vz>0) { m.vz=-m.vz; press_v += m.vz; }

		if (m.x-m.radius < 0 && m.vx<0) { m.vx=-m.vx; press_v += m.vx; }
		if (m.y-m.radius < 0 && m.vy<0) { m.vy=-m.vy; press_v += m.vy; }
		if (m.z-m.radius < 0 && m.vz<0) { m.vz=-m.vz; press_v += m.vz; }

		s.total_presure_energy+= 0.5*MASS*square(press_v);

		const double vel = std::sqrt(square(m.vx)+square(m.vy)+square(m.vz));

		// other particles?

		// do a knn search
		double query_pt[3] = { m.x, m.y, m.z };

		const size_t num_results = 4;
		size_t ret_index[num_results];
		double out_dist_sqr[num_results];
		nanoflann::KNNResultSet<double> resultSet(num_results);
		resultSet.init(ret_index, out_dist_sqr );
		index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams());

		for (size_t j=0;j<num_results;j++)
		{
			if (ret_index[j]==i)
				continue;
			if (out_dist_sqr[j]> square( 2 * RADIUS ) )
				break;

			// Compute how much they "overlap":
			const double COL_OVERLAP =  2 * RADIUS - std::sqrt( out_dist_sqr[j] );

			// Ellastic collision between particles "i1" <-> "i2"
			const size_t i1 = i;
			const size_t i2 = ret_index[j];

			if (already_collide[i2] || already_collide[i])
				continue;

			already_collide[i] = true;
			already_collide[i2] = true;

			TMass &m1 = m;
			TMass &m2 = objs[i2];

			// Direction along centers: (perpendicular)  1-->2
			TPoint3D dir_cent( m2.x-m1.x, m2.y-m1.y, m2.z-m1.z );
			const double dist = dir_cent.norm();
			if (std::abs(dist)<1e-9) continue;
			dir_cent *= 1.0/dist;

			// Are they approaching each other? Skip if they're receding.
			const TPoint3D relative_vel( m2.vx-m1.vx, m2.vy-m1.vy, m2.vz-m1.vz );
			const double dv = relative_vel.x*dir_cent.x+relative_vel.y*dir_cent.y+relative_vel.z*dir_cent.z;
			if (dv>0) continue;


			CMatrixTemplateNumeric<double>  R = mrpt::math::generateAxisBaseFromDirection(dir_cent.x,dir_cent.y,dir_cent.z);

			const TPoint3D dir_tang_u( R(0,1),R(1,1),R(2,1) );
			const TPoint3D dir_tang_v( R(0,2),R(1,2),R(2,2) );

			// Get tangent components of both velocities (dot product = proyection)
			const double vt1u = dir_tang_u.x*m1.vx + dir_tang_u.y*m1.vy + dir_tang_u.z*m1.vz;
			const double vt2u = dir_tang_u.x*m2.vx + dir_tang_u.y*m2.vy + dir_tang_u.z*m2.vz;
			const double vt1v = dir_tang_v.x*m1.vx + dir_tang_v.y*m1.vy + dir_tang_v.z*m1.vz;
			const double vt2v = dir_tang_v.x*m2.vx + dir_tang_v.y*m2.vy + dir_tang_v.z*m2.vz;

			const double vp1 = dir_cent.x*m1.vx + dir_cent.y*m1.vy + dir_cent.z*m1.vz;
			const double vp2 = dir_cent.x*m2.vx + dir_cent.y*m2.vy + dir_cent.z*m2.vz;

			// Solutions for 2 particules of equal masses:
			const double vp1_new = vp2;
			const double vp2_new = vp1;

			// set new vels:
			m1.vx = vp1_new*dir_cent.x + vt1u*dir_tang_u.x+ vt1v*dir_tang_v.x;
			m1.vy = vp1_new*dir_cent.y + vt1u*dir_tang_u.y+ vt1v*dir_tang_v.y;
			m1.vz = vp1_new*dir_cent.z + vt1u*dir_tang_u.z+ vt1v*dir_tang_v.z;

			m2.vx = vp2_new*dir_cent.x + vt2u*dir_tang_u.x+ vt2v*dir_tang_v.x;
			m2.vy = vp2_new*dir_cent.y + vt2u*dir_tang_u.y+ vt2v*dir_tang_v.y;
			m2.vz = vp2_new*dir_cent.z + vt2u*dir_tang_u.z+ vt2v*dir_tang_v.z;

		}

	}


	// Simulate linear motion:
	// ------------------------------------
	for (size_t i=0;i<N;i++)
	{
		TMass &m = objs[i];

		m.x += m.vx * At;
		m.y += m.vy * At;
		m.z += m.vz * At;
	}
}

double chi_pdf(const double x, const double s, const int k)
{
	//return std::sqrt(2/(M_PI*pow(stats.maxwell_vel_std,DOFs)))* pow(x,DOFs-1)*exp( -0.5*square(x/stats.maxwell_vel_std) );

	const double lambda_k_2 = (k==3) ? 0.5*std::sqrt(M_PI) : 1;

	return pow(2, 1-0.5*k) * pow(x,k-1) * pow(s,-k) * exp(-0.5*square(x/s) ) / lambda_k_2;
}


void compute_stats( const TMassVector &objs, TStats &s)
{
	const size_t N=objs.size();

	s.overall_Ec = 0;

	vector<double> all_Vs;
	all_Vs.reserve(N);

	s.max_vel=0;

	for (size_t i=0;i<N;i++)
	{
		const TMass &m = objs[i];

		const double v2 = square(m.vx)+square(m.vy)+square(m.vz);
		const double v = std::sqrt(v2);

		const double Ec = 0.5 * v2; // mass = 1

		mrpt::utils::keep_max(s.max_vel, v);

		all_Vs.push_back(v);
		s.overall_Ec += Ec;
	}

	s.overall_Ec*=MASS;

	const double mean_Ec = s.overall_Ec / N;

	// Temperature:
	const double k = 1.3806503e-23;  // Boltzmann constant
	s.temperature =  mean_Ec*2.0/(k*DOFs);

	// Compute theoretical Maxwell-Boltzmann distribution of velocities:
	// a^2 is the variance of the molecules velocities in each dimension:
	const double a2 = 2 * mean_Ec / (MASS*DOFs);
	s.maxwell_vel_std = std::sqrt(a2);

	// Set fixed limit for the histogram:
	const double a = s.maxwell_vel_std;
	const double theor_vel_mean = 2*a*std::sqrt(2/M_PI);
	const double theor_vel_var  = a2*(3*M_PI-8)/M_PI;

	const double limit_vel_histo = std::max(0.01, theor_vel_mean + 3*std::sqrt(theor_vel_var) );

	// Compute histogram:
	s.histogram = mrpt::math::histogram(all_Vs,0,limit_vel_histo, nBins, true /* Normalized, such that the histogram can be interpreted as a PDF */, &s.hist_centers );

}

void draw_stats( opengl::COpenGLViewportPtr & gl_extra, const TStats &stats )
{
	const size_t W = DEFAULT_WIN_WIDTH;
	const size_t H = 150;

	mrpt::utils::CImage img(W,H, CH_RGB);
	img.filledRectangle(0,0,W-1,H-1, TColor(220,220,220) );

	// Draw histogram bars:
	const size_t N = stats.histogram.size();
	const int BORDER = 15; // px
	const double BAR_W = (W - 2*BORDER)/N;

	const double max_hist = mrpt::math::maximum(stats.histogram);

	const double kx = BAR_W;
	const double ky = (H-2*BORDER) / max_hist;

	const int BAR_W_2 = (BAR_W/2)-1;

	const TColor colBars(240,20,20);
	const TColor colTheory(0,10,250);

	for (size_t i=0;i<N;i++)
	{
		const double val = stats.histogram[i];

		const int cx  = BORDER + (i+0.5)*kx;
		const int cy0 = H-BORDER;
		const int cy1 = H-BORDER - val * ky;

		img.filledRectangle( cx-BAR_W_2, cy1, cx+BAR_W_2, cy0, colBars);
	}

	int prev_cx=0,prev_cy=0;
	for (size_t i=0;i<N;i++)
	{
		const double x   = stats.hist_centers[i];
		const double val = chi_pdf(x,stats.maxwell_vel_std, DOFs);

		const int cx  = BORDER + (i+0.5)*kx;
		const int cy  = H-BORDER - val * ky;

		if (i>0)
			img.line(prev_cx,prev_cy, cx,cy, colTheory );

		prev_cx=cx;
		prev_cy=cy;
	}

	// Draw labels:
	img.selectTextFont("6x13");
	for (size_t i=0;i<N;i+=4)
	{
		const double x = stats.hist_centers[i];

		const int cx  = BORDER + (i+0.5)*kx;
		const int cy  = H-BORDER+2;

		img.textOut(cx,cy,mrpt::format("%.02f",x), TColor::black);
	}
	img.textOut(W-35,H-BORDER+2,"(m/s)", TColor::black);

	// Legend:
	img.selectTextFont("6x13");
	img.textOut(W-140,5,mrpt::format("Temperature: %.04fK",stats.temperature), TColor::black);
	img.textOut(W-200,20,"Red: Real velocities histogram", TColor::black);
	img.textOut(W-200,35,"Blue: Maxwell-Boltzmann pdf", TColor::black);


	gl_extra->setImageView_fast(img);
}


// ------------------------------------------------------
//						MAIN
// ------------------------------------------------------
int main()
{
	try
	{
		BoltzmannDemo();
		return 0;
	} catch (std::exception &e)
	{
		std::cout << "MRPT exception caught: " << e.what() << std::endl;
		return -1;
	}
	catch (...)
	{
		printf("Untyped exception!!");
		return -1;
	}
}

