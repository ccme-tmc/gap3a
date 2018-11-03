        SUBROUTINE outbwin
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%                                                                 %%
C %% This subroutine creates the output file case.oubwin             %% 
C %% which contains all the informations for the charge density      %%
C %% self-consistency.                                               %%
C %%                                                                 %%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C Definition of the variables :
C ----------------------------
        USE almblm_data
        USE common_data
        USE file_names
        USE prnt
        IMPLICIT NONE
        INTEGER :: is, ik, ou
C
        WRITE(buf,'(a)')'Writing the file case.outbwin...'
        CALL printout(0)
C
        DO is=1,ns
C ====================================
C Definition of the file case.oubwin :
C ====================================
C If the computations is spin-polarized, the output file is divided 
C in two files : case.oubwinup and case.oubwindn
          IF(ifSP.AND.is==1) THEN 
            ou=oubwinup
          ELSEIF(ifSP.AND.is==2) THEN
            ou=oubwindn
          ELSE
            ou=oubwin
          ENDIF
C =======================================
C General informations about the system :
C =======================================
C
C Number of k-points in the I-BZ
          WRITE(ou,'(i6)') nk
C Definition of the Spin-orbit flag ifSO
          IF(ifSO) THEN
           WRITE(ou,'(i6)') 1
          ELSE
           WRITE(ou,'(i6)') 0 
          ENDIF
C ====================================================
C Description of the main properties of each k-point :
C ====================================================
          DO ik=1,nk
C Description of the if-included flag
            IF(kp(ik,is)%included) THEN
             WRITE(ou,'(i6)') 1
            ELSE
             WRITE(ou,'(i6)') 0
            ENDIF
            IF(kp(ik,is)%included) THEN
C Range of bands included at each k-point
             WRITE(ou,'(2(i6))') kp(ik,is)%nb_bot,kp(ik,is)%nb_top
C Weight associated to each k-point (for the simple point integration)
             WRITE(ou,*) kp(ik,is)%weight
            ENDIF
          ENDDO   ! End of the ik loop
        ENDDO     ! End of the is loop
C
        RETURN
        END 
        
        

