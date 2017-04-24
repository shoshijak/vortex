/*
 *  tree.h
 *  ex5 solution
 *
 *  Created by Dmitry Alexeev on May 30, 2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#ifndef RUN_WITH_OMP
#include <hpx/include/lcos.hpp>
#endif


#pragma once

struct Node
{
	int part_start, part_end;
	int child_id;
	double mass, xcom, ycom, r;
	int level, morton_id;
	int parent;

        void rootsetup(int part_start, int part_end, int level, int morton_id)
	{
		this->part_start = part_start;
		this->part_end = part_end;
		this->child_id = 0;
		this->level = level;
		this->morton_id = morton_id;
	}

        void setup(int part_start, int part_end, int level, int morton_id,
                   int nodeid)
        {
                this->part_start = part_start;
                this->part_end = part_end;
                this->child_id = 0;
                this->level = level;
                this->morton_id = morton_id;
                this->parent = nodeid;
        }

};

#ifdef RUN_WITH_OMP
void build(const double* const x, const double* const y, const double* mass, const int n, const int k,
						double* xsorted, double* ysorted, double* mass_sorted,
						Node *nodes, double *expansions,
						double &exTm, double &morTm,
						double &sorTm, double &reoTm, double &bldTm, int &nnodes);
#else
hpx::future<void> build(const double* const x, const double* const y, const double* mass, const int n, const int k,
        double* xsorted, double* ysorted, double* mass_sorted,
           Node *nodes, double *expansions,
           double &exTm, double &morTm,
           double &sorTm, double &reoTm, double &bldTm, int &nnodes, std::uint64_t hpx_task_threshold);
#endif

