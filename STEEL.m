function Stress=STEEL(Strain,Properties)

%Input
%Strain: Strain value
%Properties: A vector defining steel properties as shown:
       %(1):Yield plateau (0 if there is no yield plateau)
       %(2):Ratio of stress at proportional limit to 0.2% stress
       %(3):0.2% stress
       %(4):Ultimate stress
       %(5):Ultimate strain
       %(6):Modulus of elasticity
       
%Output
%Stress: Stress value

      LAMDAT=Properties(1);    
      T3=Properties(2);     
      SIGM02=Properties(3); 
      SIGMAU=Properties(4); 
      EPSULT=Properties(5); 
      ESTEEL=Properties(6); 
      
      if LAMDAT>0             
      RNU0=SIGM02/(SIGM02+LAMDAT*ESTEEL);
      RNUULT=SIGMAU/EPSULT/ESTEEL;
      ATTA02=0.2*SIGM02/(SIGMAU-SIGM02);
      RNU02=1.2*SIGM02/(ESTEEL*(1.2*SIGM02/ESTEEL+1.5*LAMDAT));
      EPSEL=(LAMDAT+SIGM02/ESTEEL)*RNU0;
      end
      
      if LAMDAT>0  
      if abs(Strain) <= LAMDAT+SIGM02/ESTEEL
      RNU0=1.0;
      SIGMAU=SIGM02;
      SIGM02=0.99*SIGM02;
      EPSULT=LAMDAT+SIGMAU/ESTEEL;
      EPSEL=0.97*SIGMAU/ESTEEL;
      SIGMEL=0.97*SIGMAU;
      ATTA02=(SIGM02-SIGMEL)/(SIGMAU-SIGMEL);
      RNU02=SIGM02/(SIGM02+0.002*ESTEEL);
      RNUULT=SIGMAU/ESTEEL/EPSULT;
      end
      end
      
      if LAMDAT==0
      SIGMEL=T3*SIGM02;
      EPSEL=SIGMEL/ESTEEL;
      RNU0=1;
      ATTA02=(SIGM02-SIGMEL)/(SIGMAU-SIGMEL);
      RNU02=SIGM02/(SIGM02+0.002*ESTEEL);
      RNUULT=SIGMAU/ESTEEL/EPSULT;
      end
      
      E1=(RNU0-RNUULT)^2*(ATTA02^2-1)+(RNUULT-RNU02)^2;
      E1=E1/((ATTA02^2-ATTA02)*(RNU0-RNUULT)^2);
      if E1>2 
          E1=2;
      end
      E2=1-E1;
      Q9=RNUULT*EPSULT-EPSEL;
      Q1=Strain^2*E2*(RNU0-RNUULT)^2+Q9^2;
      RNU1=1;
      if (abs(Strain) >= EPSEL)
      Q2=abs(Strain)*(RNU0-RNUULT)^2*(E1*Q9-2.*EPSEL*E2)-2*RNUULT*Q9^2;
      Q3=Q9^2*(2.*RNUULT-RNU0)*RNU0-(RNU0-RNUULT)^2*EPSEL*(E1*Q9-EPSEL*E2);
      RNU1=(-Q2+sqrt(Q2^2-4*Q1*Q3))/2/Q1;
      end
      v=abs(RNU1);
      
      Stress=v*ESTEEL*Strain;
       
end
      