#ifndef hgctileboard_h
#define hgctileboard_h

#include <map>
#include <string>

namespace HGCTileBoards{

  struct TileBoard{
    // int layer;
    // int iring;
    // int iphi;
    // int toa_hits;
    // int tdc_hits;
    // int adc_hits;
    int plane; 
    std::string itype;
    int nROCs;
    int irmin;
    int irmax;
    int ir;
    bool isEngine;
  };

  const std::map<std::string, TileBoard> tb_map = {
      {"L34_J8",     {34,"J8",1,18,25,1,false}},
      {"L34_K4",     {34,"K4",1,26,29,2,true}},
      {"L35_J8",     {35,"J8",1,18,25,3,false}},
      {"L35_K6",     {35,"K6",1,26,31,4,true}},
      {"L36_J8",     {36,"J8",1,18,25,1,false}},
      {"L36_K7",     {36,"K7",1,26,33,2,true}},
      {"L37_J8",     {37,"J8",1,18,25,3,false}},
      {"L37_K7",     {37,"K7",1,26,33,4,true}},
      {"L38_C5",     {38,"C5",1,13,17,1,false}},
      {"L38_D8",     {38,"D8",1,18,25,2,false}},
      {"L38_E8",     {38,"E8",1,26,33,3,false}},
      {"L38_G3",     {38,"G3",1,34,37,4,true}},
      {"L39_C5",     {39,"C5",1,13,17,1,false}},
      {"L39_D8",     {39,"D8",1,18,25,2,false}},
      {"L39_E8",     {39,"E8",1,26,33,3,false}},
      {"L39_G4",     {39,"G4",1,34,39,4,true}},
      {"L40_B11B12", {40,"B11B12",2,6,17,1,false}},
      {"L40_D8",     {40,"D8",1,18,25,2,false}},
      {"L40_E8",     {40,"E8",1,26,33,3,false}},
      {"L40_G8",     {40,"G8",1,34,41,4,true}},
      {"L41_B11B12", {41,"B11B12",2,6,17,1,false}},
      {"L41_D8",     {41,"D8",1,18,25,2,false}},
      {"L41_E8",     {41,"E8",1,26,33,3,false}},
      {"L41_G8",     {41,"G8",1,34,41,4,true}},
      {"L42_B11B12", {42,"B11B12",2,6,17,1,false}},
      {"L42_D8",     {42,"D8",1,18,25,2,false}},
      {"L42_E8",     {42,"E8",1,26,33,3,false}},
      {"L42_G8",     {42,"G8",1,34,41,4,true}},
      {"L43_B11B12", {43,"B11B12",2,6,17,1,false}},
      {"L43_D8",     {43,"D8",1,18,25,2,false}},
      {"L43_E8",     {43,"E8",1,26,33,3,false}},
      {"L43_G8",     {43,"G8",1,34,41,4,true}},
      {"L44_A5A6",   {44,"A5A6",1,0,5,1,false}},
      {"L44_B12",    {44,"B12",2,6,17,2,false}},
      {"L44_D8",     {44,"D8",1,18,25,3,false}},
      {"L44_E8",     {44,"E8",1,26,33,4,false}},
      {"L44_G8",     {44,"G8",1,34,41,5,true}},
      {"L45_A5A6",   {45,"A5A6",1,0,5,1,false}},
      {"L45_B12",    {45,"B12",2,6,17,2,false}},
      {"L45_D8",     {45,"D8",1,18,25,3,false}},
      {"L45_E8",     {45,"E8",1,26,33,4,false}},
      {"L45_G8",     {45,"G8",1,34,41,5,true}},
      {"L46_A5A6",   {46,"A5A6",1,0,5,1,false}},
      {"L46_B12",    {46,"B12",2,6,17,2,false}},
      {"L46_D8",     {46,"D8",1,18,25,3,false}},
      {"L46_E8",     {46,"E8",1,26,33,4,false}},
      {"L46_G8",     {46,"G8",1,34,41,5,true}},
      {"L47_A5A6",   {47,"A5A6",1,0,5,1,false}},
      {"L47_B12",    {47,"B12",2,6,17,2,false}},
      {"L47_D8",     {47,"D8",1,18,25,3,false}},
      {"L47_E8",     {47,"E8",1,26,33,4,false}},
      {"L47_G6",     {47,"G6",1,34,39,5,true}},
  };
};


#endif

