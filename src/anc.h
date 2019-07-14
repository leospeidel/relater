#ifndef ANC_HPP
#define ANC_HPP

////////////////////////////////////
// Classes for Trees and AncesTrees

#include <iostream>
#include <iomanip>
#include <list>
#include <limits>
#include <deque>
#include <ctime>
#include <tgmath.h>

//#include "gzstream.h"

//Data structure of a node in a tree
struct Node{

  Node* parent      = NULL;
  Node* child_left  = NULL;
  Node* child_right = NULL;
  int label;

  float num_events     = 0.0; //on branch on top of parent
  int SNP_begin        = 0;
  int SNP_end          = 0;
  double branch_length = 0.0; //of branch on top of parent

  void operator=(const Node& node){
    parent        = node.parent;
    child_left    = node.child_left;
    child_right   = node.child_right;
    label         = node.label;
    num_events    = node.num_events;
    SNP_begin     = node.SNP_begin;
    SNP_end       = node.SNP_end;
    branch_length = node.branch_length;
  }
  bool operator==(Node node) const {
    if(node.label==label){
      return true;
    }else{
      return false;
    }
  }
  bool operator!=(Node node) const {
    if(node.label!=label){
      return true;
    }else{
      return false;
    }
  }

};


#endif //ANC_HPP
