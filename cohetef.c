#include <stdio.h>
#include <math.h>
//#include "gsl_rng.h" //Libreria para generación de números aleatorios
//Me hubiese gustado programarlo para cualquier direccion del meteorito, pero no he tenido tiempo

#define N 4 //numero de ecuaciones diferenciales
#define h 0.1 //paso fijo
#define TMAX 1000000


#define G 6.67E-11
#define Mt 5.9736E24
#define Ml 0.07349E24
#define dtl 3.844E8
#define Rt 6.37816E6
#define Rl 1.7374E6
#define w 2.6617E-6
#define PI  3.14159265359

#define RadMet 10000 
#define densidad 208   //  Kg/M^3
#define Bomba 2E19

//gsl_rng *tau;

double delta, mu;
double deltaM, muM;

double f0(double y[N]); //calcula dr/dt
double f1(double y[N]); //calcula dphy/dt
double f2(double y[N], double t); //calcula dpr/dt
double f3(double y[N], double t); //calcula dpphy/dt

void RungeKutta(double variable[N], double t);


double f0M(double y[N]); //calcula dr/dt
double f1M(double y[N]); //calcula dphy/dt
double f2M(double y[N], double t); //calcula dpr/dt
double f3M(double y[N], double t); //calcula dpphy/dt

void RungeKuttaMet(double variable[N], double t);


double dist(double x[N], double y[N]);
double min(double x, double y);

int main()
{
    //extern gsl_rng *tau;    //Puntero al estado del número aleatorio
    //int semilla=27544;     //Semilla del generador de números aleatorios

	extern double delta, mu;
    double v, theta, v0; //Valores iniciales
    double c[N]; //cohete[N]
	double t,texplosion; //corredor de tiempo e instante de la explosion
	int i, j;

	FILE *fresult;

    double M[N];    //Meteorito[N]
    double mMet;    //Masa Meteorito    
    double radioMet;    //Radio Meteorito
    double aux;

    double  tapogeo, tfrenar, tvrad;    //Instantes donde se activarán los propulsores

    double m1[N], m2[N];    // meterotito1[N], meteorito2[N]
    double L[N], T[N]; //Luna[N], Tierra[N]

    FILE *fmeteorito;

    FILE *fimpulsos;

	//y[0]=r
	//y[1]=phy
	//y[2]=p_r
	//y[3]=p_phy
       
    //tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    //gsl_rng_set(tau,semilla); //Inicializamos la semilla


	fresult=fopen("cohetef.txt", "w");  //Guardaremos posiciones Tierra, Luna, cohete
    fmeteorito=fopen("mets.txt", "w");  //Guardaremos posiciones del Met y los mets
    fimpulsos=fopen("impulsos.txt", "w");   //Guardaremos los valores de los impulsos suministrados
	

    mMet=densidad*4.*PI*RadMet*RadMet*RadMet/3.;   //Masa del meteorito
    radioMet=RadMet/dtl;    // radio del meteorito (reescalado)

    
    delta=G*Mt/(dtl*dtl*dtl);
    mu=Ml/Mt;

    deltaM=G*Mt/(dtl*dtl*dtl);
    muM=mMet/Mt;


/********************************************************************************************/
    //CONDICIONES INICIALES DE CADA OBJETO

    //Luna
    L[0]=1.;
    L[1]=0.;    

    //Tierra
    T[0]=0.;
    T[1]=0.;

    //Cohete
	c[0]=Rt*1./dtl; //r inicial (reescalada)
    c[1]=10*PI/180.; //phy inicial (latitud del lanzamiento)

	v0=1.25*sqrt(2.*G*Mt/Rt); //vel inicial 1.25*(velocidad de escape de la Tierra)   
	v=v0*1./dtl;// v inicial (reescalada)

	theta=(10-2.8)*PI/180.; //Direccion de lanzamiento del cohete (prueba y error)
    
    c[2]=v*cos(theta-c[1]);// pr inicial (reescalada, no tiene en cuanta la masa cohete)
    c[3]=c[0]*v*sin(theta-c[1]);// p_phy inicial (reescalada)

    //Meteorito
    M[0]=10.;    //Distancia inicial (reescalada)
    M[1]=-25.*PI/180.;  //Angulo de incidencia, puede ser cualquier numero [0,2*PI]

    v0=5000; //Velocidad de acercamiento del meteorito, 5000 m/s
    v=v0*1./dtl;// v inicial (reescalada)

    M[2]=-v;// pr inicial (reescalada)
    M[3]=0.;// p_phy inicial (reescalada)

/*****************************************************************************************/


    //Inicializamos así estas variables para obtener los tiempos donde activar los propulsores    
    aux=dtl;
    texplosion=TMAX;
    tfrenar=TMAX;

	j=250; //Para que guarde las condiciones iniciales
	
    //Hacemos avanzar el tiempo
    for(t=0;t<TMAX;t=t+h)
    {
		//escribo coordenadas al fichero si j=250
		if(j==250)
		{
            fprintf(fresult,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,T[0],T[0],c[0]*cos(c[1]),
c[0]*sin(c[1]),cos(L[1]),sin(L[1])); //Tierra, cohete, Luna
           
            //Si la cabeza no ha explotado guarda las posiciones de los met como la del Met grande
            if(t<=texplosion){
                fprintf(fmeteorito, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t,M[0]*cos(M[1]),M[0]*sin(M[1]),M[0]*cos(M[1]),M[0]*sin(M[1]),M[0]*cos(M[1]),M[0]*sin(M[1])); 
                }
            
            //Una vez explota estudiamos las dos partes generadas
            else{
                fprintf(fmeteorito, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t,M[0]*cos(M[1]),M[0]*sin(M[1]),m1[0]*cos(m1[1]),m1[0]*sin(m1[1]),m2[0]*cos(m2[1]),m2[0]*sin(m2[1]));
                }                
            
            j=0;
		}
		j++;		
        
        //Impulso de los motores para acelerar la nave cerca de la Luna
        if((t<=42822.91)&&(t>=42822.82))    
        {   
            fprintf(fimpulsos,"Pre impulso: c[2]\tc[3]: %lf\t%lf\n", c[2],c[3]);
            c[3]=c[3]-15.*0.0000005;
            c[2]=c[2]+6.*0.0000005;
            fprintf(fimpulsos,"Post impulso: c[2]\tc[3]: %lf\t%lf\n", c[2],c[3]);
        }
        
        //Impulso de los motores para pararse en la trayectoria del Met
        if((t<=221166.51)&&(t>=221166.49))   
        {   
            fprintf(fimpulsos,"Pre impulso: c[2]\tc[3]: %lf\t%lf\n", c[2],c[3]);
            c[3]=0.;
            c[2]=0.;
            fprintf(fimpulsos,"Post impulso: c[2]\tc[3]: %lf\t%lf\n", c[2],c[3]);
        }

        //Impulso de los motores hacia la Tierra para tener v_relativa nula 
        if((t<=324448.21)&&(t>=324448.19))   
        {   
            fprintf(fimpulsos,"Pre impulso: c[2]\tc[3]: %lf\t%lf\n", c[2],c[3]);
            c[2]=M[2]+0.000000001;
            fprintf(fimpulsos,"Post impulso: c[2]\tc[3]: %lf\t%lf\n", c[2],c[3]);
        }
        

        //Se produce la explosion si:
        if((dist(c,M)<radioMet)&&(fabs(c[2]-M[2])<0.000000009)) 
        //La velocidad relativa no es exactamente nula, pero es menos de 0.5m/s
        {
            //Una vez impacta tenemos que
            texplosion=t;
            printf("Instante del impacto: %lf\n",texplosion);

            
            //Las condiciones iniales de las 2 partes del meteorito son
            m1[0]=M[0];
            m2[0]=M[0];

            m1[1]=M[1];
            m2[1]=M[1];

            v0=sqrt(2.*Bomba/mMet); //E/2=1/2*(1/2 mMet)*velocidad_tangencial^2
            v=v0*1./dtl; //velocidad generada por al explosion (reescalada)

            m1[2]=M[2];
            m2[2]=M[2];

            m1[3]=v;
            m2[3]=-v;

            
            //Colocamos el meteorito en las posiciones iniciales (gif más bonito)
            M[0]=10.;
            M[1]=-25.*PI/180.;
            M[2]=0.;
            M[3]=0.; 
        }


        if(t<texplosion)
        {
            /***********************************************************************************/
            //Así evaluamos los distintos tiempos obtenidos para el uso de los propulsores
            //Para hacer uso de ellos hay que comentar partes de arriba, es necesario ir en 
                //progreso siguiendo la evolución temporal del sistema


            //Para vez cuando usar los propulsores para salir de la órbita lunar   
            if(dist(c,L)<aux)
            {
                aux=dist(c,L);
                tapogeo=t;
            }
            
            //Para saber cuando usar los propulsores para mantenernos en la trayectoria del
                //meteorito
            if((c[1]<=M[1]+0.0000001)&&(c[1]>=M[1]-0.0000001))
            {
                if(tfrenar>t)
                {
                    tfrenar=t;
                    printf("Instante para frenar vtang: %lf\n",tfrenar);
                }
            }
        
            //Para saber cuando usar los propulsores para tener una velocidad relativa nula con 
                //respecto al meteorito
            if(c[0]<M[0]-radioMet)
                tvrad=t;

            /*********************************************************************************/
                
            //AQUI LAS POSICIONES DE LOS OBJETOS EVOLUCIONAN TEMPORALMENTE            

            //Predomina la acción gravitatoria lunar (sobre el cohete) si:
            if(10.*10.*10.*10.*dist(c,M)>dist(c,L))
            {
                RungeKutta(c,t);
                RungeKutta(M,t);
                L[1]=w*h+L[1];  //evolución del ángulo de la Luna
            }
        
            //Predomina la acción gravitatoria del Meteorito (sobre el cohete) si:
            else
            {
                RungeKuttaMet(c,t); //El cohete es afectado por la gravedad del Meteorito
                RungeKutta(M,t);    //El Met sigue siendo afectado solo por la Tierra y Luna
                L[1]=w*h+L[1]; 
            }
        }


        if(t>=texplosion) //pongo esto en vez de else porque queda más claro
        {
            RungeKutta(m1,t);
            RungeKutta(m2,t);
            L[1]=w*h+L[1]; 
            
            //Estudiaremos si las partes impactan contra la Luna o la Tierra.
            if(dist(m1,T)<Rt/dtl)
                printf("PELIGRO, met1 IMPACTO CONTRA LA TIERRA, dist(m1,T)=%lf\n", dist(m1,T));
            if(dist(m2,T)<Rt/dtl)
                printf("PELIGRO, met2 IMPACTO CONTRA LA TIERRA, dist(m2,T)=%lf\n", dist(m2,T));
            if(dist(m1,T)<Rl/dtl)
                printf("PELIGRO, met1 IMPACTO CONTRA LA LUNA, dist(m1,L)=%lf\n", dist(m1,L));
            if(dist(m1,T)<Rl/dtl)
                printf("PELIGRO, met2 IMPACTO CONTRA LA LUNA, dist(m2,L)=%lf\n", dist(m2,L));
            
        }
    }    

    printf("Instante de distancia minima: %lf\n", tapogeo);
    printf("Instante para acelerar hacia la Tierra: %lf\n",tvrad);
    
	fclose(fresult);
    fclose(fmeteorito);
    fclose(fimpulsos);

    return 0;
}


//función que te devuelve el valor de dr/dt 
double f0(double y[N])
{
	double r;

	r=y[2];	
	return r;
}


//función que te devuelve el valor de dphi/dt 
double f1(double y[N])
{
	double phy;

	phy=y[3]/y[0]/y[0];	
	return phy;
}

//función que te devuelve el valor de dp_r/dt
double f2(double y[N], double t)
{
	extern double delta, mu;
	double pr, r2;
	
	r2=sqrt(1+y[0]*y[0]-2*y[0]*cos(y[1]-w*t));
	pr=y[3]*y[3]/(y[0]*y[0]*y[0])-delta*(1./(y[0]*y[0])+mu*(y[0]-cos(y[1]-w*t))/(r2*r2*r2));
	
	return pr;
}

//función que te devuelve el valor de dp_phy/dt
double f3(double y[N], double t)
{
	extern double delta, mu;
	double pphy, r2;
	
	r2=sqrt(1+y[0]*y[0]-2*y[0]*cos(y[1]-w*t));
	pphy=-mu*delta*y[0]*sin(y[1]-w*t)/(r2*r2*r2);
	return pphy;
}



void RungeKutta(double variable[N], double t)
{
    double  k1[N], k2[N], k3[N], k4[N];
    double aux[N];
    int i;

	//cálculo de k1
	k1[0]=h*f0(variable);
	k1[1]=h*f1(variable);
	k1[2]=h*f2(variable,t);
	k1[3]=h*f3(variable,t);

	//calculo el vector auxiliar para evaluar k2
	for(i=0;i<N;i++)
	{
		aux[i]=variable[i]+k1[i]/2.;
	}

	//cálculo de k2
	k2[0]=h*f0(aux);
    k2[1]=h*f1(aux);
	k2[2]=h*f2(aux,t+h/2.);
	k2[3]=h*f3(aux,t+h/2.);

	//calculo el vector auxiliar para evaluar k3
	for(i=0;i<N;i++)
	{
		aux[i]=variable[i]+k2[i]/2.;
	}

	//cálculo de k3
	k3[0]=h*f0(aux);
    k3[1]=h*f1(aux);
	k3[2]=h*f2(aux,t+h/2.);
	k3[3]=h*f3(aux,t+h/2.);
		
	//calculo el vector auxiliar para evaluar k4
	for(i=0;i<N;i++)
	{
		aux[i]=variable[i]+k3[i];
	}

	//cálculo de k4
	k4[0]=h*f0(aux);
    k4[1]=h*f1(aux);
	k4[2]=h*f2(aux,t+h);
	k4[3]=h*f3(aux,t+h);
		
	//Nuevas coordenadas
	for (i=0;i<N;i++)
	{
		variable[i]=variable[i]+(k1[i]+2.*k2[i]+2.*k3[i]+k4[i])/6.;	
	}

    return;
}


void RungeKuttaMet(double variable[N], double t)
{
    double  k1[N], k2[N], k3[N], k4[N];
    double aux[N];
    int i;

	//cálculo de k1
	k1[0]=h*f0M(variable);
	k1[1]=h*f1M(variable);
	k1[2]=h*f2M(variable,t);
	k1[3]=h*f3M(variable,t);

	//calculo el vector auxiliar para evaluar k2
	for(i=0;i<N;i++)
	{
		aux[i]=variable[i]+k1[i]/2.;
	}

	//cálculo de k2
	k2[0]=h*f0M(aux);
    k2[1]=h*f1M(aux);
	k2[2]=h*f2M(aux,t+h/2.);
	k2[3]=h*f3M(aux,t+h/2.);

	//calculo el vector auxiliar para evaluar k3
	for(i=0;i<N;i++)
	{
		aux[i]=variable[i]+k2[i]/2.;
	}

	//cálculo de k3
	k3[0]=h*f0M(aux);
    k3[1]=h*f1M(aux);
	k3[2]=h*f2M(aux,t+h/2.);
	k3[3]=h*f3M(aux,t+h/2.);
		
	//calculo el vector auxiliar para evaluar k4
	for(i=0;i<N;i++)
	{
		aux[i]=variable[i]+k3[i];
	}

	//cálculo de k4
	k4[0]=h*f0M(aux);
    k4[1]=h*f1M(aux);
	k4[2]=h*f2M(aux,t+h);
	k4[3]=h*f3M(aux,t+h);
		
	//Nuevas coordenadas
	for (i=0;i<N;i++)
	{
		variable[i]=variable[i]+(k1[i]+2.*k2[i]+2.*k3[i]+k4[i])/6.;	
	}

    return;
}

//función que te devuelve el valor de dr/dt 
double f0M(double y[N])
{
	double r;

	r=y[2];	
	return r;
}


//función que te devuelve el valor de dphi/dt 
double f1M(double y[N])
{
	double phy;

	phy=y[3]/y[0]/y[0];	
	return phy;
}

//función que te devuelve el valor de dp_r/dt
double f2M(double y[N], double t)
{
	extern double deltaM, muM;
	double pr, r2;
	
	r2=sqrt(1+y[0]*y[0]-2*y[0]*cos(y[1]-w*t));
	pr=y[3]*y[3]/(y[0]*y[0]*y[0])-deltaM*(1./(y[0]*y[0])+muM*(y[0]-cos(y[1]-w*t))/(r2*r2*r2));
	
	return pr;
}

//función que te devuelve el valor de dp_phy/dt
double f3M(double y[N], double t)
{
	extern double deltaM, muM;
	double pphy, r2;
	
	r2=sqrt(1+y[0]*y[0]-2*y[0]*cos(y[1]-w*t));
	pphy=-muM*deltaM*y[0]*sin(y[1]-w*t)/(r2*r2*r2);
	return pphy;
}


double dist(double x[N], double y[N])  //Calcula la distancia entre 2 cuerpos
{ 
    double d1,d2;
    double norma;

    // Las componentes del vector distancia serán
    d1=x[0]*cos(x[1])-y[0]*cos(y[1]);
    d2=x[0]*sin(x[1])-y[0]*sin(y[1]);

    norma=sqrt(d1*d1+d2*d2);

    return norma;
}

double min(double x, double y)
{
    double minimo;

    if(x<y)
        minimo=x;

    else
        minimo=y;

    return minimo;
}

