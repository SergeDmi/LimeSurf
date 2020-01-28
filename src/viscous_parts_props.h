/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   elastic_parts_prop.h
 * Author: dmitrief
 *
 * Created on July 9, 2018, 11:16 AM
 */
#include <math.h>
#include <sstream>
//#include "glossary.h"
#include "Aboria.h"
#include "spin_parts_props.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

using namespace Aboria;
#ifndef VISCOUS_PARTS_PROPS_H
#define VISCOUS_PARTS_PROPS_H


class Viscous_set_props : public Part_set_props 
{
    friend class Viscous_part_set;
    
    public:
        /// constructor
        Viscous_set_props();
        //Viscous_set_props(const Glossary );
        Viscous_set_props(const YAML::const_iterator );
        void init();
        void Read_config(const YAML::const_iterator );
    protected:
        // prestrain : how much the imported link is stretched
        //double prestrain;
        //int power_law;
        double R0;              // Equilibrium distance between parts
        double minR,minR2;      // Minimum interaction distance (avoid dvg)
        double Rmax;            // Max interaction distance
        double Rsearch;         // Search radius
        //double pressure;        // A pressure
        double k_rep,p_rep;     // repulsive potential and power law
        double k_att,p_att;     // attractive potential and power law
        double k_bend,p_bend;   // bending potential and power law
        double p_align; // alignment potential and power law
        double viscosity;
        double L;               // box size
        
        double renorm_rate;     // rate at which normals are normalized can be quite low
        double Fmax,Fmax2;      // Maximum force
        double relax;           // linkers relaxation rate
        vdouble3 corner_0,corner_1; // Corners of box
        int max_neighbours;
        
        //double x_max, y_max, z_max ;
        //double x_conf,y_conf,z_conf;
        
        // @todo ; with keylist
        int init_shape;
        int init_number;
        double init_radius;
        
        //bool elastic;
        //double k_elast;
    

};


#endif /* VISCOUS_PARTS_PROP_H */
