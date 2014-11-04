C***********************************************************************
C
C     diagonalizzazione matrice reale simmetrica: libreria EISPACK
C
      subroutine eispd(md,ndim,a,e,c,w)
C
C     MD    = max dim
C     NDIM  = dimension
C     A     = matrix to be diag.
C     E     = eigenvalues
C     C     = eigenvectors; A and C may coincide. If distinct, A is
C             unaltered
C     W     = work array
C
      integer md,ndim,ierr
      double precision a(md,*),c(md,*),w(md,*),e(*)
      ierr=0
      call tred2(md,ndim,a,e,w,c)
      call tql2(md,ndim,e,w,c,ierr)
      if (ierr.ne.0) then 
         write(6,*) ' EISPACK: diag. failed;   IERR =',ierr
         stop 12
      end if
      return
      end
C***********************************************************************
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                            
C                                                                               
      INTEGER I,J,K,L,N,II,NM,JP1                                               
C     REAL*8 A(NM,N),D(N),E(N),Z(NM,N)                                          
      REAL*8 A(NM,*),D(*),E(*),Z(NM,*)                                          
      REAL*8 F,G,H,HH,SCALE                                                     
      REAL*8 DSQRT,DABS,DSIGN                                                   
C                                                                               
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,            
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.           
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).           
C                                                                               
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A                      
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING                       
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                                    
C                                                                               
C     ON INPUT:                                                                 
C                                                                               
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL                 
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM                  
C          DIMENSION STATEMENT;                                                 
C                                                                               
C        N IS THE ORDER OF THE MATRIX;                                          
C                                                                               
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE                  
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.                       
C                                                                               
C     ON OUTPUT:                                                                
C                                                                               
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;            
C                                                                               
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL                 
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;              
C                                                                               
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                        
C          PRODUCED IN THE REDUCTION;                                           
C                                                                               
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.                    
C                                                                               
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,                
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY                 
C                                                                               
C     ------------------------------------------------------------------        
C                                                                               
      DO 100 I = 1, N                                                           
C                                                                               
         DO 100 J = 1, I                                                        
            Z(I,J) = A(I,J)                                                     
  100 CONTINUE                                                                  
C                                                                               
      IF (N .EQ. 1) GO TO 320                                                   
C     :::::::::: FOR I=N STEP -1 UNTIL 2 DO -- ::::::::::                       
      DO 300 II = 2, N                                                          
         I = N + 2 - II                                                         
         L = I - 1                                                              
         H = 0.0D0                                                              
         SCALE = 0.0D0                                                          
         IF (L .LT. 2) GO TO 130                                                
C     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::               
         DO 120 K = 1, L                                                        
  120    SCALE = SCALE + DABS(Z(I,K))                                           
C                                                                               
         IF (SCALE .NE. 0.0D0) GO TO 140                                        
  130    E(I) = Z(I,L)                                                          
         GO TO 290                                                              
C                                                                               
  140    DO 150 K = 1, L                                                        
            Z(I,K) = Z(I,K) / SCALE                                             
            H = H + Z(I,K) * Z(I,K)                                             
  150    CONTINUE                                                               
C                                                                               
         F = Z(I,L)                                                             
         G = -DSIGN(DSQRT(H),F)                                                 
         E(I) = SCALE * G                                                       
         H = H - F * G                                                          
         Z(I,L) = F - G                                                         
         F = 0.0D0                                                              
C                                                                               
         DO 240 J = 1, L                                                        
            Z(J,I) = Z(I,J) / H                                                 
            G = 0.0D0                                                           
C     :::::::::: FORM ELEMENT OF A*U ::::::::::                                 
            DO 180 K = 1, J                                                     
  180       G = G + Z(J,K) * Z(I,K)                                             
C                                                                               
            JP1 = J + 1                                                         
            IF (L .LT. JP1) GO TO 220                                           
C                                                                               
            DO 200 K = JP1, L                                                   
  200       G = G + Z(K,J) * Z(I,K)                                             
C     :::::::::: FORM ELEMENT OF P ::::::::::                                   
  220       E(J) = G / H                                                        
            F = F + E(J) * Z(I,J)                                               
  240    CONTINUE                                                               
C                                                                               
         HH = F / (H + H)                                                       
C     :::::::::: FORM REDUCED A ::::::::::                                      
         DO 260 J = 1, L                                                        
            F = Z(I,J)                                                          
            G = E(J) - HH * F                                                   
            E(J) = G                                                            
C                                                                               
            DO 260 K = 1, J                                                     
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)                          
  260    CONTINUE                                                               
C                                                                               
  290    D(I) = H                                                               
  300 CONTINUE                                                                  
C                                                                               
  320 D(1) = 0.0D0                                                              
      E(1) = 0.0D0                                                              
C     :::::::::: ACCUMULATION OF TRANSFORMATION MATRICES ::::::::::             
      DO 500 I = 1, N                                                           
         L = I - 1                                                              
         IF (D(I) .EQ. 0.0D0) GO TO 380                                         
C                                                                               
         DO 360 J = 1, L                                                        
            G = 0.0D0                                                           
C                                                                               
            DO 340 K = 1, L                                                     
  340       G = G + Z(I,K) * Z(K,J)                                             
C                                                                               
            DO 360 K = 1, L                                                     
               Z(K,J) = Z(K,J) - G * Z(K,I)                                     
  360    CONTINUE                                                               
C                                                                               
  380    D(I) = Z(I,I)                                                          
         Z(I,I) = 1.0D0                                                         
         IF (L .LT. 1) GO TO 500                                                
C                                                                               
         DO 400 J = 1, L                                                        
            Z(I,J) = 0.0D0                                                      
            Z(J,I) = 0.0D0                                                      
  400    CONTINUE                                                               
C                                                                               
  500 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
C     :::::::::: LAST CARD OF TRED2 ::::::::::                                  
      END                                                                       
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                          
C                                                                               
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR                                     
C     REAL*8 D(N),E(N),Z(NM,N)                                                  
      REAL*8 D(*),E(*),Z(NM,*)                                                  
      REAL*8 B,C,F,G,H,P,R,S,MACHEP                                             
      REAL*8 DSQRT,DABS,DSIGN                                                   
C                                                                               
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,             
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND             
C     WILKINSON.                                                                
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).           
C                                                                               
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS                    
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.                       
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO                      
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                          
C     FULL MATRIX TO TRIDIAGONAL FORM.                                          
C                                                                               
C     ON INPUT:                                                                 
C                                                                               
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL                 
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM                  
C          DIMENSION STATEMENT;                                                 
C                                                                               
C        N IS THE ORDER OF THE MATRIX;                                          
C                                                                               
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;                  
C                                                                               
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX                
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;                       
C                                                                               
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE                   
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS              
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN                
C          THE IDENTITY MATRIX.                                                 
C                                                                               
C      ON OUTPUT:                                                               
C                                                                               
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN                  
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT                  
C          UNORDERED FOR INDICES 1,2,...,IERR-1;                                
C                                                                               
C        E HAS BEEN DESTROYED;                                                  
C                                                                               
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC                   
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,             
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED               
C          EIGENVALUES;                                                         
C                                                                               
C        IERR IS SET TO                                                         
C          ZERO       FOR NORMAL RETURN,                                        
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN                       
C                     DETERMINED AFTER 30 ITERATIONS.                           
C                                                                               
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,                
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY                 
C                                                                               
C     ------------------------------------------------------------------        
C                                                                               
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING             
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.           
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC                
C                ON S360 ::::::::::                                             
      parameter (machep=16.d0**(-13))
c     DATA MACHEP/Z'3410000000000000'/
C                                                                               
      IERR = 0                                                                  
      IF (N .EQ. 1) GO TO 1001                                                  
C                                                                               
      DO 100 I = 2, N                                                           
  100 E(I-1) = E(I)                                                             
C                                                                               
      F = 0.0D0                                                                 
      B = 0.0D0                                                                 
      E(N) = 0.0D0                                                              
C                                                                               
      DO 240 L = 1, N                                                           
         J = 0                                                                  
         H = MACHEP * (DABS(D(L)) + DABS(E(L)))                                 
         IF (B .LT. H) B = H                                                    
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::                 
         DO 110 M = L, N                                                        
            IF (DABS(E(M)) .LE. B) GO TO 120                                    
C     :::::::::: E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT                       
C                THROUGH THE BOTTOM OF THE LOOP ::::::::::                      
  110    CONTINUE                                                               
C                                                                               
  120    IF (M .EQ. L) GO TO 220                                                
  130    IF (J .EQ. 30) GO TO 1000                                              
         J = J + 1                                                              
C     :::::::::: FORM SHIFT ::::::::::                                          
         L1 = L + 1                                                             
         G = D(L)                                                               
         P = (D(L1) - G) / (2.0D0 * E(L))                                       
         R = DSQRT(P*P+1.0D0)                                                   
         D(L) = E(L) / (P + DSIGN(R,P))                                         
         H = G - D(L)                                                           
C                                                                               
         DO 140 I = L1, N                                                       
  140    D(I) = D(I) - H                                                        
C                                                                               
         F = F + H                                                              
C     :::::::::: QL TRANSFORMATION ::::::::::                                   
         P = D(M)                                                               
         C = 1.0D0                                                              
         S = 0.0D0                                                              
         MML = M - L                                                            
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::                     
         DO 200 II = 1, MML                                                     
            I = M - II                                                          
            G = C * E(I)                                                        
            H = C * P                                                           
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150                              
            C = E(I) / P                                                        
            R = DSQRT(C*C+1.0D0)                                                
            E(I+1) = S * P * R                                                  
            S = C / R                                                           
            C = 1.0D0 / R                                                       
            GO TO 160                                                           
  150       C = P / E(I)                                                        
            R = DSQRT(C*C+1.0D0)                                                
            E(I+1) = S * E(I) * R                                               
            S = 1.0D0 / R                                                       
            C = C * S                                                           
  160       P = C * D(I) - S * G                                                
            D(I+1) = H + S * (C * G + S * D(I))                                 
C     :::::::::: FORM VECTOR ::::::::::                                         
            DO 180 K = 1, N                                                     
               H = Z(K,I+1)                                                     
               Z(K,I+1) = S * Z(K,I) + C * H                                    
               Z(K,I) = C * Z(K,I) - S * H                                      
  180       CONTINUE                                                            
C                                                                               
  200    CONTINUE                                                               
C                                                                               
         E(L) = S * P                                                           
         D(L) = C * P                                                           
         IF (DABS(E(L)) .GT. B) GO TO 130                                       
  220    D(L) = D(L) + F                                                        
  240 CONTINUE                                                                  
C     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::                  
      DO 300 II = 2, N                                                          
         I = II - 1                                                             
         K = I                                                                  
         P = D(I)                                                               
C                                                                               
         DO 260 J = II, N                                                       
            IF (D(J) .GE. P) GO TO 260                                          
            K = J                                                               
            P = D(J)                                                            
  260    CONTINUE                                                               
C                                                                               
         IF (K .EQ. I) GO TO 300                                                
         D(K) = D(I)                                                            
         D(I) = P                                                               
C                                                                               
         DO 280 J = 1, N                                                        
            P = Z(J,I)                                                          
            Z(J,I) = Z(J,K)                                                     
            Z(J,K) = P                                                          
  280    CONTINUE                                                               
C                                                                               
  300 CONTINUE                                                                  
C                                                                               
      GO TO 1001                                                                
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN                              
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::                      
 1000 IERR = L                                                                  
 1001 RETURN                                                                    
C     :::::::::: LAST CARD OF TQL2 ::::::::::                                   
      END                                                                       
