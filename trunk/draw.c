/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "draw.h"


/*********************************************************/

void Draw_Tree(edge *b_root, tdraw *w, arbre *tree)
{
  Init_Tdraw_Struct(w);
  Get_Tree_Box_Width(w,tree);
  Dist_To_Root(b_root,tree);
  w->max_dist_to_root = Get_Max_Dist_To_Root(tree);
  Get_X_Coord(b_root,w,tree);
  Get_Y_Coord(b_root,w,tree);
}

/*********************************************************/

void Print_Postscript_Header(int n_pages, FILE *fp)
{
  fprintf(fp,"%%!PS-Adobe-3.0\n");
  fprintf(fp,"%%%%DocumentFonts: Times-Roman Times-Roman\n");
  fprintf(fp,"%%%%Creator: Stephane Guindon\n");
  fprintf(fp,"%%%%Title: tree\n");
  fprintf(fp,"%%%%BeginFeature: *PageSize\n"); 
  fprintf(fp,"a4\n");
  fprintf(fp,"%%%%EndFeature\n");
  fprintf(fp,"%%%%EndComments\n");
  fprintf(fp,"%%%%Pages: %d\n",n_pages);

  fprintf(fp,"/lt {lineto} bind def\n");
  fprintf(fp,"/mt {moveto} bind def\n");
  fprintf(fp,"/sc {setrgbcolor} bind def\n");

  fprintf(fp,"/clipbox\n");
  fprintf(fp,"{\n");
  fprintf(fp,"newpath\n");
  fprintf(fp,"40 40 moveto\n");
  fprintf(fp,"560 40 lineto\n");
  fprintf(fp,"560 810 lineto\n");
  fprintf(fp,"40 810 lineto\n");
  fprintf(fp,"40 40 lineto\n");
  fprintf(fp,"closepath\n");
  fprintf(fp,"clip\n");
  fprintf(fp,"} bind def\n");
  
  fprintf(fp,"/Times-Roman findfont\n");
  fprintf(fp,"12 scalefont\n");
  fprintf(fp,"setfont\n");


}

/*********************************************************/

void Print_Postscript_EOF(FILE *fp)
{
  fprintf(fp,"%%%%Trailer\n");
  fprintf(fp,"%%%%EOF\n");
}

/*********************************************************/

void Print_Tree_Postscript(edge *b_root, FILE *fp, int tree_num, tdraw *w, arbre *tree)
{
  int i;
  

  fprintf(fp,"%%%%Page: %d %d\n",tree_num+1,tree_num+1); 
  fprintf(fp,"clipbox\n");
  fprintf(fp,"stroke\n");
  fprintf(fp,"50 50 translate\n");
  fprintf(fp,"newpath\n");

  printf("\n. Reading branch %d, has probability %f attached to it.",b_root->num,b_root->prob_sel_regime);

  if(b_root->prob_sel_regime <= 0.1)
    fprintf(fp,".0 .0 1. sc\n");
  else if(b_root->prob_sel_regime > 0.1 && b_root->prob_sel_regime <= 0.2)
    fprintf(fp,".0 .5 1. sc\n");
  else if(b_root->prob_sel_regime > 0.2 && b_root->prob_sel_regime <= 0.3)
    fprintf(fp,".0 1. 1. sc\n");
  else if(b_root->prob_sel_regime > 0.3 && b_root->prob_sel_regime <= 0.4)
    fprintf(fp,".0 1. .5 sc\n");
  else if(b_root->prob_sel_regime > 0.4 && b_root->prob_sel_regime <= 0.5)
    fprintf(fp,".0 1. .0 sc\n");
  else if(b_root->prob_sel_regime > 0.5 && b_root->prob_sel_regime <= 0.6)
    fprintf(fp,".5 1. .0 sc\n");
  else if(b_root->prob_sel_regime > 0.6 && b_root->prob_sel_regime <= 0.7)
    fprintf(fp,"1. 1. 0. sc\n");
  else if(b_root->prob_sel_regime > 0.7 && b_root->prob_sel_regime <= 0.8)
    fprintf(fp,"1. .5 0. sc\n");
  else if(b_root->prob_sel_regime > 0.8 && b_root->prob_sel_regime <= 0.9)
    fprintf(fp,"1. 0. 0. sc\n");
  else if(b_root->prob_sel_regime > 0.9)
    fprintf(fp,"1. .0 .0 sc\n");


  fprintf(fp,"%d %d mt\n",w->xcoord[b_root->left->num],w->ycoord[b_root->left->num]);
  fprintf(fp,"%d %d lt\n",0,w->ycoord[b_root->left->num]);
  fprintf(fp,"%d %d lt\n",0,w->ycoord[b_root->rght->num]);
  fprintf(fp,"%d %d lt\n",w->xcoord[b_root->rght->num],w->ycoord[b_root->rght->num]);
  fprintf(fp,"stroke\n");

  fprintf(fp,"%d %d mt\n",w->xcoord[b_root->left->num],w->ycoord[b_root->left->num]);
  if(b_root->left->tax) 
    {
      fprintf(fp,"(%s) show\n",b_root->left->name);
    }
  else
    {
      For(i,3)
	if((b_root->left->v[i]) && (b_root->left->v[i] != b_root->rght))
	  Print_Tree_Postscript_Pre(b_root->left,b_root->left->v[i],fp,w,tree);
    }

  fprintf(fp,"%d %d mt\n",w->xcoord[b_root->rght->num],w->ycoord[b_root->rght->num]);

  if(b_root->rght->tax) 
    {
      fprintf(fp,"(%s) show\n",b_root->rght->name);
    }
  else
    {
      For(i,3)
	if((b_root->rght->v[i]) && (b_root->rght->v[i] != b_root->left))
	  Print_Tree_Postscript_Pre(b_root->rght,b_root->rght->v[i],fp,w,tree);
    }

  fprintf(fp,"closepath\n");
  fprintf(fp,"stroke\n");
  fprintf(fp,"showpage\n");
}

/*********************************************************/

void Print_Tree_Postscript_Pre(node *a, node *d, FILE *fp, tdraw *w, arbre *tree)
{
  int i;
  double min,max;
  double step;

  min = 0.0;
  max = 4.;

  step = (max-min)/13.;

  fprintf(fp,"gsave\n");
  
  For(i,3)
    if(a->v[i] == d)
      {
	printf("\n. Reading branch %d, has probability %f attached to it.",a->b[i]->num,a->b[i]->prob_sel_regime);
       
	 if(a->b[i]->prob_sel_regime <= min+1.*step)
	  fprintf(fp,".0 1. 1. sc\n");
	else if(a->b[i]->prob_sel_regime > min+1.*step && a->b[i]->prob_sel_regime <= min+2.*step)
	  fprintf(fp,".0 1. .8 sc\n");
	else if(a->b[i]->prob_sel_regime > min+2.*step && a->b[i]->prob_sel_regime <= min+3.*step)
	  fprintf(fp,".0 1. .5 sc\n");
	else if(a->b[i]->prob_sel_regime > min+3.*step && a->b[i]->prob_sel_regime <= min+4.*step)
	  fprintf(fp,".0 1. .3 sc\n");
	else if(a->b[i]->prob_sel_regime > min+4.*step && a->b[i]->prob_sel_regime <= min+5.*step)
	  fprintf(fp,".0 1. .0 sc\n");
	else if(a->b[i]->prob_sel_regime > min+5.*step && a->b[i]->prob_sel_regime <= min+6.*step)
	  fprintf(fp,".25 1. 0. sc\n");
	else if(a->b[i]->prob_sel_regime > min+6.*step && a->b[i]->prob_sel_regime <= min+7.*step)
	  fprintf(fp,".5 1. 0. sc\n");
	else if(a->b[i]->prob_sel_regime > min+7.*step && a->b[i]->prob_sel_regime <= min+8.*step)
	  fprintf(fp,".75 1. .0 sc\n");
	else if(a->b[i]->prob_sel_regime > min+8.*step && a->b[i]->prob_sel_regime <= min+9.*step)
	  fprintf(fp,"1. 1. .0 sc\n");
	else if(a->b[i]->prob_sel_regime > min+9.*step && a->b[i]->prob_sel_regime <= min+10.*step)
	  fprintf(fp,"1. .75 .0 sc\n");
	else if(a->b[i]->prob_sel_regime > min+10.*step && a->b[i]->prob_sel_regime <= min+11.*step)
	  fprintf(fp,"1. .5 .0 sc\n");
	else if(a->b[i]->prob_sel_regime > min+11.*step && a->b[i]->prob_sel_regime <= min+12.*step)
	  fprintf(fp,"1. .25 .0 sc\n");
	else if(a->b[i]->prob_sel_regime > min+12.*step && a->b[i]->prob_sel_regime <= min+13.*step)
	  fprintf(fp,"1. .0 0. sc\n");
	break;
      }

  fprintf(fp,"%d %d mt\n",w->xcoord[a->num],w->ycoord[a->num]);
  fprintf(fp,"%d %d lt\n",w->xcoord[a->num],w->ycoord[d->num]);
  fprintf(fp,"%d %d lt\n",w->xcoord[d->num],w->ycoord[d->num]);

  if(d->tax) 
    {
      fprintf(fp,"(%s) show \n",d->name);
      fprintf(fp,"stroke\n");
      fprintf(fp,"grestore\n");
      return;
    }
  else
    {
      fprintf(fp,"stroke\n");
      fprintf(fp,"grestore\n");
      For(i,3)
	if(d->v[i] != a) Print_Tree_Postscript_Pre(d,d->v[i],fp,w,tree);
    }


  return;
}

/*********************************************************/

void Dist_To_Root_Pre(node *a, node *d, edge *b, arbre *tree)
{
  int i;

  d->dist_to_root = a->dist_to_root + b->l;

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a) Dist_To_Root_Pre(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/

void Dist_To_Root(edge *b_root, arbre *tree)
{  
  int i;

  b_root->left->dist_to_root = b_root->l / 2.;
  b_root->rght->dist_to_root = b_root->l / 2.;

  For(i,3)
    {
      if((b_root->left->v[i]) && (b_root->left->v[i] != b_root->rght))
	Dist_To_Root_Pre(b_root->left,b_root->left->v[i],b_root,tree);
      
      if((b_root->rght->v[i]) && (b_root->rght->v[i] != b_root->left))
	Dist_To_Root_Pre(b_root->rght,b_root->rght->v[i],b_root,tree);
    }
}

/*********************************************************/

void Get_X_Coord_Pre(node *a, node *d, edge *b, tdraw *w, arbre *tree)
{
  int i;

  w->xcoord[d->num] =  d->dist_to_root * (double)w->tree_box_width/w->max_dist_to_root;

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a) Get_X_Coord_Pre(d,d->v[i],d->b[i],w,tree);
    }
}

/*********************************************************/

void Get_X_Coord(edge *b_root, tdraw *w, arbre *tree)
{
  int i;

  w->xcoord[b_root->left->num] = b_root->left->dist_to_root * (double)w->tree_box_width/w->max_dist_to_root;
  w->xcoord[b_root->rght->num] = b_root->rght->dist_to_root * (double)w->tree_box_width/w->max_dist_to_root;

  For(i,3)
    {
      if((b_root->left->v[i]) && (b_root->left->v[i] != b_root->rght))
	Get_X_Coord_Pre(b_root->left,b_root->left->v[i],b_root, w, tree);

      if((b_root->rght->v[i]) && (b_root->rght->v[i] != b_root->left))
	Get_X_Coord_Pre(b_root->rght,b_root->rght->v[i],b_root, w, tree);
    }
}

/*********************************************************/

void Get_Y_Coord(edge *b_root, tdraw *w, arbre *tree)
{
  int next_y_slot;
  next_y_slot = 0;
  Get_Y_Coord_Post(b_root->left,b_root->rght,b_root,&next_y_slot,w,tree);
  Get_Y_Coord_Post(b_root->rght,b_root->left,b_root,&next_y_slot,w,tree);
}

/*********************************************************/

void Get_Y_Coord_Post(node *a, node *d, edge *b, int *next_y_slot, tdraw *w, arbre *tree)
{
  int i;

  if(d->tax) 
    {
      w->ycoord[d->num] = *next_y_slot + (int)(w->page_height / (2.*tree->n_otu));
      (*next_y_slot) += (int)(w->page_height / (tree->n_otu));
    }
  else
    {
      int d1, d2;

      d1 = d2 = -1;
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Get_Y_Coord_Post(d,d->v[i],d->b[i],next_y_slot,w,tree);
	      if(d1<0) d1 = i;
	      else     d2 = i;
	    }
	}
      w->ycoord[d->num] = (w->ycoord[d->v[d1]->num] + w->ycoord[d->v[d2]->num])/2.; 
    }
}

/*********************************************************/

tdraw *Make_Tdraw_Struct(arbre *tree)
{
  tdraw *w;

  w = (tdraw *)mCalloc(1,sizeof(tdraw));
  w->xcoord = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  w->ycoord = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));

  return w;
}

/*********************************************************/

void Init_Tdraw_Struct(tdraw *w)
{
  w->page_width  = 510;
  w->page_height = 770;
}

/*********************************************************/

void Get_Tree_Box_Width(tdraw *w, arbre *tree)
{
  int i;
  int max_name_len, curr_len;

  max_name_len = curr_len = 0;
  For(i,tree->n_otu)
    {
      curr_len = (int)strlen(tree->noeud[i]->name);
      if(curr_len > max_name_len) max_name_len = curr_len;
    }

  w->tree_box_width = w->page_width - max_name_len * 8.66667;
}

/*********************************************************/

double Get_Max_Dist_To_Root(arbre *tree)
{
  double mx;
  int i;

  mx = .0;
  For(i,tree->n_otu)
    {
      if(tree->noeud[i]->dist_to_root > mx)
	{
	  mx = tree->noeud[i]->dist_to_root;
	}
    }

  return mx;
}

/*********************************************************/
/*********************************************************/
