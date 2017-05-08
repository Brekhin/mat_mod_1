 
        WRITE (*,'(A)')'Input end interval'
        READ *, B
        WRITE(*,'(A)')'input step'
        READ *, H

        OPEN(1, FILE='1')

        XI = 0
        YI_P1 = 0
        YI_P2 = 0
        YI_P3 = 0
        YI_P4 = 0
        YI_RK = 0
        POQ = 0
        YI_EXPL = 0

        DO WHILE (XI.LE.(B - H/2))
            YI_RK = RUNGE_KUTTA(XI, YI_RK, H)   
            POQ = POQED(XI, POQ, H)
            YI_EXPL = EXPL_SCHEME(XI, YI_EXPL, H)
            
            XI = XI + H

            YI_P1 = PIKAR_1(XI)
            YI_P2 = PIKAR_2(XI)
            YI_P3 = PIKAR_3(XI)
            YI_P4 = PIKAR_4(XI)
            
            WRITE (1, 100) XI, YI_P1, YI_P2, YI_P3, YI_P4, YI_RK, POQ, YI_EXPL
            100 FORMAT (8f10.4)
        END DO 
    
        CLOSE(1)
 
        PAUSE
        STOP
    END

    FUNCTION PIKAR_1(X)
        PIKAR_1 = X**3/3
    END
    
    FUNCTION PIKAR_2(X)
        PIKAR_2 = X**3/3*(1 + X**4/21)
    END
    
    FUNCTION PIKAR_3(X)
        PIKAR_3 = X**3 / 3*(1 + X**4 * (1 / 21 + X**4 * (2 / 693 + X**4 / 19845)))
    END
    
    FUNCTION PIKAR_4(X)
        A = X**4
        B = A * (4 / 1113959385 + A / 36625634325.0)
        C = A * (82 / 12442815 + A * (662/3479404005.0 + B))
        D = A * (13 / 72765 + C)
        E = A * (2 / 693 + D)
        F = A * (1 / 21 + E)
        G = X**3 / 3
        PIKAR_4 = G * (1 + F)
    END
    
    FUNCTION RUNGE_KUTTA(XI, YI, H)
        HALF = H / 2
        YI_HALF = HALF * (XI**2 + YI**2) + YI
        YI_DHALF = (XI + HALF)**2 + YI_HALF**2
        RUNGE_KUTTA = YI + H * YI_DHALF
    END
    
    
     FUNCTION POQED(XI, POQ, H)
        B = H * XI + H**2
        A = B**2 + H * POQ 
        C = 1 - 4 * A
        IF (C >= 0) THEN
            POQA = (1 - sqrt(C)) / (H * 2)
        ELSE
            POQA = -1
        END IF
        POQED = POQA    
    END
    

    FUNCTION EXPL_SCHEME(XI, YI, H)
        EXPL_SCHEME = YI + H * ((XI**2) + (YI**2))
    END

