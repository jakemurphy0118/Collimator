#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "Visualize.hh"

using namespace CLHEP;

Visualize::Visualize()
{}

Visualize::~Visualize()
{}

void Visualize::CoordinateAxes()
{
  //Redraw the coordinate axis
  G4Polyline x_axis;
  //Set red line color                                                      
  G4Colour red(1.0,0.0,0.0);
  G4VisAttributes att(red);
  x_axis.push_back( G4Point3D (0.,0.,0.));
  x_axis.push_back( G4Point3D (5.*cm,0.,0.));
  x_axis.SetVisAttributes(att);


}
