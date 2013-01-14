CC = gcc
LIBS = -lm
# LIBS = -L /usr/local/lib/ -lm -lglut -lGLU -lGL -lIL -lILU -lILUT 
# LIBS = -L /usr/local/lib/ -lm -lglut -lGLU -lGL
# PROG = POSSELSITEID
# EXEC = pos
PROG = FITMODEL
EXEC = fitmodel
# PROG = TREEVIEW
# EXEC = treeview
# PROG = EVOLVE
# EXEC = evolve
# PROG = KARIN
# EXEC = karin
# PROG = RF
# EXEC = rf

DFLAG =  -DUNIX -D$(PROG) 
CFLAGS = -O2 -fomit-frame-pointer -Wall $(DFLAG)
# CFLAGS = -g -Wall $(DFLAG)
# CFLAGS = -Wall $(DFLAG) 
# OBJS = main.o utilities.o optimiz.o lk.o models.o free.o options.o eigen.o eigenmb.o gl2ps.o
OBJS = main.o utilities.o optimiz.o lk.o models.o free.o options.o eigen.o eigenmb.o draw.o

$(EXEC) : $(OBJS)
	$(CC) -o $(EXEC) $(OBJS) $(LIBS)

clean :
	@rm -f *.o 




######################################################################################################
