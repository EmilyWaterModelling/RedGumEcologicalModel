!**************************************************************************************************************
!
!                          River Red Gum Ecological Response Model (RRGERM)
!
!  RRGERM version 01 developed by E.J. Barbour October 2015 with contributions in study design from 
!  G. Kuczera, B. Croke, and P.D Driver.

!  Results from this code are provided in the manuscript (under review): 
!  "Uncertainty in wetland ecological response for environmental flow management"
!
!*************************************************************************************************************



MODULE processMod

IMPLICIT NONE
SAVE
  INTEGER :: i,k,confail,totnumirr, numrows
  REAL::sumEcoScore,sumEcoScoreup,sumEcoScorelow, summhealth, sumallirr,avEcoScoreup,avEcoScorelow,avEcoScore,minEcoScore

  TYPE dateType
   INTEGER :: day, month, year
  END TYPE dateType

  TYPE objectiveType
   CHARACTER(10) :: dateStr
   TYPE (dateType) :: date
   REAL :: Div1, Div3, flow,Bool, ecoScore,ecoScoreup,ecoScorelow,ratefall,WyVol,GWDepth,GWDepth2, rain,ecoScoreupExp(6), ecoScorelowExp(6)
  END TYPE objectiveType
  TYPE (objectiveType), ALLOCATABLE :: obj(:)

  TYPE GCSeventType
    INTEGER::classif,Dur, TotDur
  END TYPE GCSeventType
  TYPE (GCSeventType), ALLOCATABLE::GCSevent(:)


! Check if in optimisation mode 
! When not in optimisation, write out intermediate values for testing and reporting purposes
! 1= optimisation, 2= write out results ****************************
  INTEGER, PARAMETER :: opt=1

  CONTAINS


SUBROUTINE getFlowRainValues
  IMPLICIT NONE
  CHARACTER:: objective(50), header(100)
  CHARACTER(20) :: dummy, fmt
  CHARACTER(80) :: flowfile, rainfile
  REAL :: P1
  INTEGER :: j


  IF (opt==1) THEN   !In optimisation mode - used modelled flow at Booligal from IQQM
    flowfile = 'ConMod\BoolFlow.txt'
    rainfile = 'ConMod\75007_06Rainfall_1898_2006.csv'
    numrows = 39658
  ELSE               !Not in optimisation mode - use observed flow at Booligal
    flowfile = 'ConMod\BooligalObsGapsFilled_Jul1953_Jun2013.csv'
    rainfile = 'ConMod\OxleyRainObsGapsFilled_Jul1953_Jun2013.csv'
    numrows = 21915
  END IF

  OPEN (unit=70, file=flowfile, status='old')
  OPEN (unit=71,file=rainfile, status='old')
  
  ALLOCATE (obj(numrows))
  ALLOCATE (GCSevent(numrows))
   
  DO i = 1, 3
     READ(70,*)
     READ(71,*)
  END DO

  DO k=1,numrows
     fmt = '(10x,f12.0)'      
     READ(70,'(a,1x,f12.0)') obj(k)%dateStr, obj(k)%Bool  
     READ(71,fmt) obj(k)%rain

    !This was a leftover from previous versions - at some point change all to one or the other for consistency.
     obj(k)%flow=obj(k)%Bool
  END DO 
  END SUBROUTINE getFlowRainValues
  
  
  SUBROUTINE estGWDepth
    IMPLICIT NONE
    REAL, PARAMETER:: tau=1950, m=0.007500, c=17.85700,tau2=1950,m2=0.0037,c2=8.9286
    REAL :: a,a2
    INTEGER::j
    REAL(KIND=8):: Qinit,Qinit2
    REAL(KIND=8) :: Qadj(numrows),Qadj2(numrows)
  !------------------------------------------------------------------------------
  !This subroutine estimates groundwater depth based on a regression relationship between flow at Booligal
  !  and groundwater levels at gauge GW036721 in the Great Cumbung Swamp. The regression relationship involves
  !  a Nash cascade of two storages: Qadj(i)=2a*Qadj(i-1) - a^2*Qadj(i-2) + (1-a)^2*Q(t), a=e^(-1/tau)
  !  Qadj(t) = adjusted flow at time (i)
  !  tau     = time constant = 1950 (calibrated parameter)
  !  i       = time step

  !The adjusted flow Qadj(t) is then plotted against observed groundwater levels, and a linear regression relationship
  !  is derived where GWDepth=m*Qadj+c (m=slope, c=y intercept)
  !------------------------------------------------------------------------------
  

    a=exp(-1/tau)
    i=1
    DO i=1,numrows
      IF (i==1) THEN
        Qinit=(c-12.67000)/(a*m)
        Qadj(i)=a*Qinit+(1-a)*obj(i)%Bool !Set initial GW level to the average of observed values (-12.67) by setting initial adjusted flow to (c-12.67)/a*m
      ELSE IF (i==2) THEN
        Qadj(i)=a*Qadj(i-1)+(1-a)*obj(i)%Bool
      ELSE  
        Qadj(i)=2*a*Qadj(i-1)-(a**2)*Qadj(i-2)+((1-a)**2)*obj(i)%Bool   !Calculate adjusted flow time series at Booligal using an auto-regressive model Qadj(i)=2a*Qadj(i-1) - a^2*Qadj(i-2) + (1-a)^2*Q(t), a=e^(-1/tau)
      END IF
      obj(i)%GWDepth=m*Qadj(i)-c
      WRITE(40,*) obj(i)%dateStr, obj(i)%GWDepth
    END DO
    
    a2=exp(-1/tau2)
    j=1
    DO j=1,numrows
      IF (j==1) THEN
        Qinit2=(c2-6.34)/(a2*m2)
        Qadj2(j)=a2*Qinit2+(1-a2)*obj(j)%Bool !Set initial GW level to the average of observed values (-12.67) by setting initial adjusted flow to (c-12.67)/a*m
      ELSE IF (j==2) THEN
        Qadj2(j)=a2*Qadj(j-1)+(1-a2)*obj(j)%Bool      
      ELSE
        Qadj2(j)=2*a2*Qadj(j-1)-(a2**2)*Qadj(j-2)+((1-a2)**2)*obj(j)%Bool
      END IF
      obj(j)%GWDepth2=m2*Qadj2(j)-c2
    
    END DO

    !Print outputs if not in optimisation mode
    
    !GWDepth is estimated based on a linear equation between observed GW data at gauge GW036721 and Qadj
    !GWDepth2 halves the observed GW levels at gauge GW036721 and then fits a linear equation with Qadj.
    IF(opt==2) THEN
      OPEN(unit=40,file='ConMod\GWDepthEst.txt', status='unknown')
      OPEN(unit=42,file='ConMod\GWDepthEst2.txt', status='unknown')
      WRITE(40,'(a,1x,f12.4)')  (obj(i)%dateStr, obj(i)%GWDepth, i=1,numrows)
      WRITE(42,'(a,1x,f12.4)') (obj(j)%dateStr, obj(j)%GWDepth2, j=1,numrows)
      CLOSE (unit=40)
      CLOSE (unit=42)
    END IF
  END SUBROUTINE estGWDepth
  
   
 !----------------------------------------------------------------------------------
 !INUNDATION MODEL
 !This part of the code estimates the inundation duration of the GCS.
 !This is done by considering rainfall ponding (which only operates if flow is below thresholds)
 !and flow inundation based on flow at Booligal exceeding magnitude and duration thresholds
 !Firstly, wet and dry periods are identified based on flow and rainfall.
 !Secondly, these period are grouped into total wet and dry period for the GCS (given inundation factor may
 !cause wet periods to be followed by wet periods).


 !User: Set 0 to ignore rainfall, and 1 to include rainfall
 !----------------------------------------------------------------------------------
  SUBROUTINE InundationMod (maxk, GWCode, GWDepthCode, x0, RRGArea, RRGFlag)
  IMPLICIT NONE
    INTEGER ::RainOn, RainCount, FlowInundCount
    INTEGER:: GWCode, GWDepthCode
    INTEGER :: Store, StoreRRG, DryDur, StWetDur, Dur, AllRRGDryDur, InundDur, WetEnd, EventDur
    INTEGER :: j, RRGArea, RRGFlag
    INTEGER :: maxj, maxk, SumGCSevent
    REAL:: RainThresh, EffRain, InitialLoss, ContLoss, InfilParam, Infilt, PondingDepth
    REAL :: FlowThresh, FlowDurThresh, MaintFlowThresh, factor
    REAL :: x0
    REAL ::DroughtBreak, InundThresh, AllRRGInundThresh
   
    TYPE eventType
      INTEGER::classif,Dur,TotDur,InundDur,WetEnd
    END TYPE eventType
    TYPE (eventType), ALLOCATABLE::event(:)

    ALLOCATE (event(numrows))
   ! 

    !Contains user inputs for hydrology parameters
    OPEN (unit=880,file='ConMod\HydrologyParamFile.txt', status='old')
   
    !Open files only if not in optimisation mode
    IF (opt==2) THEN
      OPEN (unit=400,file='ConMod\Events.txt',status='unknown')
      OPEN (unit=500,file='ConMod\GCSEvents.txt',status='unknown')
    END IF

    ! Read in parameters defining inundation duration based on rainfall, flow and groundwater
   READ(880,*) FlowThresh      !Flow magnitude at which the GCS starts to inundate (ML/d)
   READ(880,*) FlowDurThresh   !Length of time required at the flow threshold for inundation to occur (days)
   READ(880,*) MaintFlowThresh !Flow magnitude at which GCS remains inundated after event started (ML/d)
   READ(880,*) factor          !Duration converstion factor: GCS inundation duration = factor x duration above threshold (excluding FlowDurThresh)
   READ(880,*) RainThresh      !Rainfall intensity required for inundation to occur (mm/d)
   READ(880,*) InitialLoss     !Rainfall initial loss applied at start of each rainfall inundation event (mm/d)
   READ(880,*) ContLoss        !Rainfall continuing loss (mm/d)
   READ(880,*) InfilParam      !Factor influcing the rate at which ponded rainfall infiltrates (from 0 to 1)
   READ(880,*) GWCode          !Code influencing whether access to groundwater is considered (0=off, 1=on)
   READ(880,*) RainOn          !Code influencing whether rainfall inundation is considered
   READ(880,*) RRGFlag         !Identifies whether flow thresholds are for the full RRG area (select 0) or not (select 1)
   READ(880,*) GWDepthCode     !Code influencing whether groundwater access is based on observed levels or halved
   READ(880,*) x0              !Inflexion point for GW relationship, related to maximum root depth
   CLOSE(unit=880)

   IF (InitialLoss>RainThresh) THEN
     OPEN (unit=881, file='ConMod\ERMWarnings.txt',status='unknown')
     WRITE (881,*) 'Warning: Initial rainfall loss', InitialLoss, 'greater than Rain Threshold', RainThresh, '. Check parameter file'
     CLOSE (unit=881)
   END IF

   !Initialise variables
   Dur=0                       !Wet period duration once flow exceeds magnitude and duration thresholds, or rainfall ponding occurs
   DryDur=0                    !Dry period duration when flow and rain is below thresholds
   stWetDur=0                  !Counter for when flow exceeds magnitude threshold, but not necessarily duration threshold. This is for flow only
   Store=0                     !This is a bucket concept for flow which exceeds magnitude threshold but not duration (days above threshold). 
                               !It increments to a maximum of the duration threshold, and during dry periods decreases 1day/time step to zero.
   
   StoreRRG=0
   RainCount=0                 !Counter to identify start of rainfall event and subtract inital loss.
   PondingDepth=0 
   FlowInundCount=0
   j=1

   DO i=1, numrows
     event(i)%InundDur=0
   END DO
     
   IF (RRGArea==1 .and. RRGFlag==1) THEN   !If calculating the inundation for the full RRG area (ie RRGArea=1) and the inundation for a smaller area (ie RRGFlag=1)
     FlowThresh=2700
     MaintFlowThresh=2700
     FlowDurThresh=30
   END IF

  !Loop through each time step (numrows) 
   DO i=1,numrows
   !Estimate ponding due to effective rainfall
   IF (RainOn==1) THEN                                  !Include rainfall inundation
     IF (obj(i)%rain>RainThresh) THEN                   !Start rainfall inundation
       RainCount=RainCount+1
       IF (RainCount==1) THEN
         EffRain=max(obj(i)%rain-InitialLoss,0.0)
       ELSE
         EffRain=max(obj(i)%rain-ContLoss,0.0)
       END IF

       Infilt= (PondingDepth + EffRain)*InfilParam      !Calculate infiltration and ponding depth
       PondingDepth=PondingDepth + EffRain-Infilt
     ELSE IF (PondingDepth>1) THEN
       EffRain=max(0.0,obj(i)%rain-ContLoss)
       Infilt= (PondingDepth+EffRain)*InfilParam
       PondingDepth=PondingDepth +EffRain-Infilt
     ELSE
       PondingDepth=0.0
     END IF
     IF (opt==2) THEN   !Not in optimisation model
       OPEN (unit=410,file='ConMod\PondingDepth.txt',status='unknown')
       WRITE(410,*)i, PondingDepth
     END IF
   END IF

   !Categorise wet/dry periods in the GCS based on flow thresholds
   !This firstly estimates the FLOW duration based on flows exceeding the magnitude and time thresholds
   !   and then uses the specified multiplication factor to estimate the INUNDATION duration (which is longer)
   !The initial If statement considers whether a rainfall event has just finished, yet flow is over the threshold in 
   !terms of magnitude but not yet for duration.

   IF ((Dur==0 .and. obj(i)%flow<FlowThresh) .OR. (obj(i)%flow<MaintFlowThresh .and. Dur>0)) THEN
     Store=max(Store-1,0)
     IF (PondingDepth<1.0) THEN   !Rainfall only considered if there is no flow - otherwise flow assumed to dominate       
         RainCount=0
         DryDur=DryDur+1     !Dry period

       IF (Dur>0) THEN     !This is the start of the dry period, set duration of FLOW event        
         event(j)%Dur=Dur
       !  IF (obj(i-1)%flow>MaintFlowThresh) THEN  !This is to distinguish between a dry period following a flow inundation and a rain inundation
         IF (FlowInundCount>0) THEN    
           InundDur=min(factor*(FlowInundCount),REAL(numrows-i))                     !If it follows a flow inundation, a factor is applied to extend the inundation duration in the GCS
           IF(Dur==FlowInundCount) THEN
             WetEnd=(i-1)+(InundDur-FlowInundCount)                           !Calculate the time step at which the inundation period ends
           ELSE
             WetEnd=(i-1)+(InundDur-FlowInundCount) !Changed Dur to FlowInundCount 26/10/2014
             InundDur=Dur-FlowInundCount+InundDur
           END IF
           IF (event(j)%InundDur>0) THEN                                             !check if a previous flow event extends beyond the most recent event
             IF (event(j)%WetEnd>WetEnd) THEN
               event(j)%InundDur=InundDur + event(j)%WetEnd-WetEnd    !Changed from event(j)%InundDur=event(j)%InundDur on 28/10/2014
               event(j)%WetEnd=event(j)%WetEnd
             ELSE
               event(j)%InundDur=InundDur
               event(j)%WetEnd=WetEnd
             END IF
           ELSE
               event(j)%InundDur=InundDur
               event(j)%WetEnd=WetEnd
           END IF                  
         ELSE
           IF (event(j)%WetEnd>i) THEN   !Flow inundation event was longer than rainfall event, durations already calculated
             event(j)%InundDur=min(event(j)%WetEnd-event(j-1)%TotDur,numrows-event(j-1)%TotDur)
             event(j)%WetEnd=event(j)%WetEnd
           ELSE            
             event(j)%InundDur=Dur                  !If it is a rainfall inundation, the duration has already been calculated and a flow conversion factor is not required.
             event(j)%WetEnd=i-1                    !Rainfall event so ends at previous time step
           END IF
         END IF

         IF (j>2) THEN
           IF (event(j-1)%InundDur==0 .and. event(j)%WetEnd<event(j-2)%WetEnd) THEN  !Check if previous wet event was longer than the current event
           event(j)%WetEnd=event(j-2)%WetEnd
           END IF
         END IF

         IF (j==1) THEN
           event(j)%TotDur=event(j)%Dur
         ELSE
           event(j)%TotDur=event(j)%Dur+event(j-1)%TotDur
         END IF
         IF (opt==2) THEN    !If not in optimisation mode
           WRITE(400,*)event(j)%classif,event(j)%Dur,event(j)%TotDur,event(j)%InundDur,event(j)%WetEnd           
         END IF
         Dur=0                                     !Re-set wet period duration and start of flow counter>threshold to zero
         stWetDur=0
         j=j+1         
       END IF    !End if Dur>0
       event(j)%classif=1      !Dry event
       stWetDur=0
       FlowInundCount=0
     ELSE   !Rain inundation event   (Else if ponding depth>=1)
       IF (FlowInundCount>0) THEN    !calculate duration of flow inundation in case this extends beyond the rainfall event
         EventDur=FlowInundCount
         InundDur=min(factor*FlowInundCount,real(numrows-i))
         WetEnd=(i-1)+(InundDur-FlowInundCount)         !Calculate the time step at which the inundation period ends
         IF (event(j)%InundDur>0) THEN                !check if a previous flow event extends beyond the most recent event
             IF (event(j)%WetEnd>WetEnd) THEN
               event(j)%InundDur=event(j)%InundDur
               event(j)%WetEnd=event(j)%WetEnd
             ELSE
               event(j)%Dur=EventDur
               event(j)%InundDur=InundDur
               event(j)%WetEnd=WetEnd
             END IF
           ELSE
               event(j)%Dur=EventDur
               event(j)%InundDur=InundDur
               event(j)%WetEnd=WetEnd
           END IF
       END IF  

       IF (DryDur>0) THEN
         event(j)%Dur=DryDur
         IF (j==1) THEN
           event(j)%TotDur=event(j)%Dur
         ELSE
           event(j)%TotDur=event(j)%Dur+event(j-1)%TotDur
         END IF
         IF (j==1) THEN
           event(j)%InundDur=event(j)%Dur
         ELSE IF (event(j)%TotDur<event(j-1)%WetEnd)THEN
           event(j)%InundDur=0
         ELSE
           event(j)%InundDur=event(j)%TotDur-event(j-1)%WetEnd
         END IF
         IF (opt==2) THEN
           WRITE(400,*)event(j)%classif,event(j)%Dur,event(j)%TotDur,event(j)%InundDur,event(j)%WetEnd
         END IF
         DryDur=0
         j=j+1         
       END IF   
       Dur=Dur+1
       event(j)%classif=2
       stWetDur=0                        !Added 20/7/2014 to re-set flow counter if only rainfall event
       FlowInundCount=0
     END IF  !(End if ponding depth<1)
   ELSE      !(Else if Dur==0...)
     !Flow inundation event
     stWetDur=stWetDur+1
     IF (stWetDur==1) THEN  !Start of wet period but not inundation.Calculate duration threshold based on time since previous wet event
       DroughtBreak=DryDur/(2.0*365.0)
       InundThresh=FlowDurThresh+DroughtBreak-Store  !Threshold for 50% of RRG to be wet
       AllRRGInundThresh=30+DroughtBreak-StoreRRG                    !Flow duration threshold for 100% RRG area to be wet
     END IF
     
     !A bucket concept is used to increment the number of days above threshold flow to the duration threshold.
     !This then decreased during dry period, but if only a short dry period occurs, the threshold duration is decreased
     !by the number of days left in the bucket.
     
     Store=min(real(Store+1),FlowDurThresh) 
     IF (obj(i-1)%flow>FlowThresh) THEN
       StoreRRG=min(real(StoreRRG+1),AllRRGInundThresh)
     END IF
                
     IF (stWetDur<InundThresh.and.DryDur>0 .and. PondingDepth<1) THEN     !Continue to stay in dry period until exceeds duration threshold
       DryDur=DryDur+1
       AllRRGDryDur=AllRRGDryDur+1
     ELSE IF (PondingDepth<1 .and. Dur>0 .and. stWetDur<InundThresh) THEN  !End of rain inundation event before start flow event
       event(j)%Dur=Dur
       event(j)%InundDur=Dur
       IF (j==1) THEN                           !Calculate the time step at which the inundation period ends
         event(j)%WetEnd=event(j)%InundDur
       ELSE IF (j==2) THEN
         event(j)%WetEnd=event(j-1)%TotDur+event(j)%InundDur
       ELSE  
         IF (event(j-1)%InundDur==0) THEN
           event(j)%WetEnd=i-1
           IF (event(j)%WetEnd<event(j-2)%WetEnd) THEN
             event(j)%WetEnd=event(j-2)%WetEnd
           END IF
         ELSE
           event(j)%WetEnd=event(j-2)%WetEnd+event(j-1)%InundDur+event(j)%InundDur  
         END IF
       END IF
       IF (j==1) THEN
         event(j)%TotDur=event(j)%Dur
       ELSE
         event(j)%TotDur=event(j)%Dur+event(j-1)%TotDur
       END IF
       IF (opt==2) THEN    !If not in optimisation mode
         WRITE(400,*)event(j)%classif,event(j)%Dur,event(j)%TotDur,event(j)%InundDur,event(j)%WetEnd           
       END IF
       Dur=0                                     !Re-set wet period duration >threshold to zero
       j=j+1         
       event(j)%classif=1      !Dry event
       DryDur=1
     ELSE
       IF (DryDur>0) THEN                                  !This is the start of an inundation period rain or flow
         event(j)%Dur=DryDur
         IF (j==1) THEN
           event(j)%TotDur=event(j)%Dur
         ELSE
           event(j)%TotDur=event(j)%Dur+event(j-1)%TotDur  !Determine the time step at which the previous dry period ended
         END IF
         IF (j==1) THEN
           event(j)%InundDur=event(j)%Dur                  !InundDur refers to the duration of the wet or dry period after multiplying flow by the inundation factor
         ELSE IF (event(j)%TotDur<event(j-1)%WetEnd)THEN   !Check if the end of the dry period occurs before the end of the prior inundation period
           event(j)%InundDur=0
         ELSE
           event(j)%InundDur=event(j)%TotDur-event(j-1)%WetEnd   !Otherwise calculate the dry period by reducing the length of time previously inundated
         END IF
         IF (opt==2) THEN
           WRITE(400,*)event(j)%classif,event(j)%Dur,event(j)%TotDur,event(j)%InundDur,event(j)%WetEnd
         END IF
         DryDur=0
         j=j+1         
       END IF
       Dur=Dur+1
       event(j)%classif=2
       IF (stWetDur>=InundThresh) THEN
         FlowInundCount=FlowInundCount+1
       ELSE
         FlowInundCount=0
       END IF
     END IF
   END IF
   END DO

   !This section of code is to process the last event, given in the code above, events are only processed once the next event is reached
  !(and for the last event there is no additional event).
  
   IF(event(j)%classif==1) THEN
     event(j)%Dur=DryDur
     event(j)%TotDur=event(j)%Dur+event(j-1)%TotDur
     IF (j==1) THEN
       event(j)%InundDur=0
     ELSE IF (event(j)%TotDur<event(j-1)%WetEnd)THEN
       event(j)%InundDur=0
     ELSE
       event(j)%InundDur=event(j)%TotDur-event(j-1)%WetEnd
     END IF
   ELSE IF (event(j)%classif==2) THEN
     event(j)%Dur=Dur
     event(j)%InundDur=Dur  !Not multiplied by factor as will otherwise extend beyond simulation duration
     IF (j==1) THEN
       event(j)%WetEnd=event(j)%InundDur
     ELSE IF (j==2) THEN
       event(j)%WetEnd=event(j-1)%TotDur+event(j)%InundDur
     ELSE  
       event(j)%WetEnd=min(event(j-2)%WetEnd+event(j-1)%InundDur+event(j)%InundDur,numrows) !If the previous dry event was zero, 
     END IF
   END IF
   IF (j==1) THEN
     event(j)%TotDur=event(j)%Dur
   ELSE
     event(j)%TotDur=event(j)%Dur+event(j-1)%TotDur
   END IF
   IF (opt==2) THEN
     WRITE(400,*)event(j)%classif,event(j)%Dur,event(j)%TotDur,event(j)%InundDur,event(j)%WetEnd
     CLOSE (unit=400)
   END IF

 
 !---------------------------------------------------------------------
 !Aggregate events identifed into total wet and dry periods for the GCS.
 !This is due to the inundation factor causing wet periods to coincide even if the flow fell below the threshold for short periods

 !---------------------------------------------------------------------

   maxj=j
   k=0
   SumGCSevent=0
  
   DO j=1,maxj
     IF(event(j)%classif==1) THEN  !dry period
       IF(event(j)%InundDur>0) THEN !dry period actually occurred
         k=k+1
         GCSevent(k)%classif=event(j)%classif
         GCSevent(k)%Dur=event(j)%InundDur
         SumGCSevent=SumGCSevent+GCSevent(k)%Dur
         IF (k==1) THEN
           GCSevent(k)%TotDur=GCSevent(k)%Dur
         ELSE
           GCSevent(k)%TotDur=GCSevent(k-1)%TotDur+GCSevent(k)%Dur
         END IF
       END IF
     ELSE      !wet period
       IF (j<3) THEN
         k=k+1
         GCSevent(k)%classif=event(j)%classif     
         GCSevent(k)%Dur=event(j)%InundDur
         SumGCSevent=SumGCSevent+GCSevent(k)%Dur
         IF (k==1) THEN
           GCSevent(k)%TotDur=GCSevent(k)%Dur
         ELSE
           GCSevent(k)%TotDur=GCSevent(k-1)%TotDur+GCSevent(k)%Dur
         END IF
       ELSE IF (event(j-1)%InundDur==0) THEN                !No dry period before 
         IF((SumGCSevent+event(j)%InundDur)>numrows) THEN   !Check if total wet period adds to more than simulation period
           GCSevent(k)%Dur=GCSevent(k)%Dur+(max(0,numrows-SumGCSevent))
           EXIT
         ELSE
           GCSevent(k)%Dur=event(j)%WetEnd-GCSevent(k-1)%TotDur
           GCSevent(k)%TotDur=GCSevent(k-1)%TotDur+GCSevent(k)%Dur
           SumGCSevent=event(j)%WetEnd   !Originally SumGCSevent=SumGCSevent+event(j)%InundDur but if inund=0 there is probably overlap between events
         END IF                                            !Hence their durations can't be added.
       ELSE
         k=k+1
         GCSevent(k)%classif=event(j)%classif 
         IF (event(j)%WetEnd<=numrows) THEN                   
           GCSevent(k)%Dur=event(j)%InundDur
           SumGCSevent=SumGCSevent+GCSevent(k)%Dur
           GCSevent(k)%TotDur=GCSevent(k-1)%TotDur+GCSevent(k)%Dur
         ELSE   ! too wet period extends beyond simulation duration
           GCSevent(k)%Dur=event(j)%InundDur - (event(j)%WetEnd-numrows)
           EXIT
         END IF
       END IF
     END IF
   END DO
 
   maxk=k
 
   IF(opt==2) THEN
     DO k=1,maxk
       WRITE(500,*)GCSevent(k)%classif,GCSevent(k)%Dur,GCSevent(k)%TotDur
     END DO
     CLOSE(unit=500)
   END IF

  END SUBROUTINE InundationMod

  !-------------------------------------------------------------------------
  !SUBROUTINE COMPARE INUNDATION AREA
  !Calculate GCS inundation for RRG area using a flow threshold of 2700ML/d and duration of 30 days
  !This is done to check whether the reed bed is actually wetter for longer periods if flows are much higher than 700ML/d

  !---------------------------------------------------------------------------
  SUBROUTINE CompareInundationArea (x0, maxk, GWCode, GWDepthCode)
  IMPLICIT NONE
    INTEGER :: w, s, m
    INTEGER :: maxk, maxm, RRGArea, RRGFlag, GWCode, GWDepthCode, GCSDryCount, GCSWetCount
    REAL :: RRGEvent(numrows),GCSWetDry(numrows), x0
    
    TYPE GCSeventRRGType
      INTEGER::classif,Dur
    END TYPE GCSeventRRGType
    TYPE (GCSeventRRGType), ALLOCATABLE::GCSeventRRG(:)

    ALLOCATE (GCSeventRRG(numrows))
 
 !Calculates inundation for full RRG area
 RRGArea=1
 CALL InundationMod(maxm, GWCode, GWDepthCode, x0, RRGArea, RRGFlag) 

 !GCSevent(maxm)%Dur=2
 IF (RRGFlag==0) THEN  !i.e. only calculating full RRG area so no more calculations needed other than checking for short events
   maxk=maxm
   m=1
   k=1
   DO s=1, maxm
     GCSevent(k)%classif=GCSevent(m)%classif
     IF (m<maxm) THEN
       IF (GCSevent(m+1)%Dur<=5) THEN   !If duration of NEXT event is less than or equal to 5, aggregate the next two events with the current
         IF (m==maxm-2) THEN
           GCSevent(k)%Dur=GCSevent(m)%Dur+GCSevent(m+1)%Dur+GCSevent(m+2)%Dur   !End of simulation
           EXIT
         ELSE IF (m==maxm-1) THEN
           GCSevent(k)%Dur=GCSevent(m)%Dur+GCSevent(m+1)%Dur  !End of simulation
           EXIT
         ELSE  !event not at end of simulation
           GCSevent(k)%Dur=GCSevent(m)%Dur+GCSevent(m+1)%Dur+GCSevent(m+2)%Dur
           k=k+1
           m=m+3
         END IF
       ELSE
         GCSevent(k)%Dur=GCSevent(m)%Dur
         k=k+1
         m=m+1
       END IF
     ELSE
       GCSevent(k)%Dur=GCSevent(m)%Dur   !maxm
       EXIT
     END IF   !end if m<maxm
   END DO

   maxk=k
   IF(opt==2) THEN
     OPEN (unit=500,file='ConMod\GCSEvents.txt',status='unknown')
     DO k=1,maxk
       WRITE(500,*)GCSevent(k)%classif,GCSevent(k)%Dur
     END DO
     CLOSE(unit=500)
   END IF
 END IF

 IF (RRGFlag==1) THEN   !Calculate for area of inundation selected by user (i.e. reed bed)
   DO m=1, maxm
     GCSeventRRG(m)%classif=GCSevent(m)%classif
     GCSeventRRG(m)%Dur= GCSevent(m)%Dur
   END DO

   ! Calculate GCS inundation for user input thresholds
   RRGArea=0
   CALL InundationMod(maxk, GWCode, GWDepthCode, x0, RRGArea, RRGFlag)

   !Take the classification of inundation events to generate a time series for each time step
   ! and a code specifying if the area is wet (code=1) or dry (code=0)
   w=1 
   DO m=1, maxm                               !Calculate wet/dry time series for full RRG areas
     IF (GCSeventRRG(m)%classif==1) THEN
       DO s=1, GCSeventRRG(m)%Dur
         RRGEvent(w)=0
         w=w+1
       END DO
     ELSE IF (GCSeventRRG(m)%classif==2) THEN
       DO s=1, GCSeventRRG(m)%Dur
         RRGEvent(w)=1
         w=w+1
       END DO
     END IF
   END DO

   w=1
   DO k=1, maxk                                      !Calculate wet/dry series for reed bed
     IF (GCSevent(k)%classif==1) THEN
       DO s=1, GCSevent(k)%Dur
         GCSWetDry(w)=0
         w=w+1
       END DO
     ELSE IF (GCSevent(k)%classif==2) THEN
       DO s=1, GCSevent(k)%Dur
         GCSWetDry(w)=1
         w=w+1
       END DO
     END IF
   END DO

   !Compare inundation duration for user input area and RRGArea (the full RRG area in the GCS), and make sure if RRGArea is wet, so is the user input area
   ! Then recombine the time series to summarise the duration of each wet and dry event for use in the ecosystem model
   k=1
   GCSDryCount=0
   GCSWetCount=0
   DO w=1, numrows
     IF (RRGEvent(w)==1) THEN
       GCSWetDry(w)=1
     END IF
     IF (GCSWetDry(w)==0) THEN    !Dry
       IF (w==1) THEN
         GCSDryCount=1
       ELSE 
         GCSDryCount=GCSDryCount+1
         IF (GCSWetDry(w)/=GCSWetDry(w-1)) THEN
           IF (GCSWetCount <= 5) THEN                               !If there is less than 5 days between wet/dry periods, aggregate into a single event
             GCSevent(k-1)%Dur=GCSevent(k-1)%Dur+GCSWetCount        !Otherwise in the ERM, excessively long wet periods won't be captured if a single dry day occurs
           ELSE IF (k>1) THEN
             IF (GCSevent(k-1)%classif==2) THEN
               GCSevent(k-1)%Dur = GCSevent(k-1)%Dur + GCSWetCount
             ELSE 
               GCSevent(k)%classif=2
               GCSevent(k)%Dur=GCSWetCount  
               k=k+1
             END IF
           ELSE
             GCSevent(k)%classif=2
             GCSevent(k)%Dur=GCSWetCount  
             k=k+1
           END IF
           GCSWetCount=0
         END IF
       END IF
     ELSE                             !Wet
       IF (w==1) THEN
         GCSWetCount=1
       ELSE 
         GCSWetCount=GCSWetCount+1
         IF (GCSWetDry(w)/=GCSWetDry(w-1)) THEN 
           IF (GCSDryCount <= 5) THEN
             GCSevent(k-1)%Dur=GCSevent(k-1)%Dur+GCSDryCount
           ELSE IF (k>1) THEN
             IF(GCSevent(k-1)%classif==1) THEN
               GCSevent(k-1)%Dur = GCSevent(k-1)%Dur + GCSDryCount
             ELSE
               GCSevent(k)%classif=1
               GCSevent(k)%Dur=GCSDryCount       
               k=k+1
             END IF
           ELSE
             GCSevent(k)%classif=1
             GCSevent(k)%Dur=GCSDryCount       
             k=k+1
           END IF
           GCSDryCount=0
         END IF
       END IF
     END IF
     IF (w==numrows) THEN              !If last time step, calculate the duration of the last event
       IF (GCSWetCount>0) THEN
         GCSevent(k)%classif=1
         GCSevent(k)%Dur=GCSWetCount
       ELSE
         GCSevent(k)%classif=1
         GCSevent(k)%Dur=GCSDryCount
       END IF
     END IF
   END DO

   maxk=k
   IF(opt==2) THEN
     OPEN (unit=500,file='ConMod\GCSEvents.txt',status='unknown')
     DO k=1,maxk
       WRITE(500,*)GCSevent(k)%classif,GCSevent(k)%Dur
     END DO
     CLOSE(unit=500)
   END IF
 END IF !End if RRGFlag=1

  END SUBROUTINE CompareInundationArea


  !---------------------------------------------------------------------------------------------
  !SUBROUTINE ENVFLOWPROCESSOR
  


  !---------------------------------------------------------------------------------------------
  
  SUBROUTINE envFlowProcessor
  USE ifPort
  IMPLICIT NONE
  INTEGER :: i, m, d, y, e, k, t, maxk
  INTEGER :: test, GWCode, GWDepthCode, GWCount, GWEndCount,GWCountAbove9, DryGW 
  INTEGER :: numexpertupper, numexpertlower, numexperts, numlines, numexperttoowet, TooWetCountUp(6), TooWetCountLow(6),DryGWExp(6),DryGWlowExp(6), WetGWExp(6)
  REAL :: PLInUp, PLConUp, PLInLow, PLConLow, testscore,DurAdj, DurAdjustWet, DurAdjustDry, GWAdj, GWav,num
  REAL:: f,c,startEcoScore,ScoreAdjUp,ScoreAdjLow,GWscale,startScoreup,startScorelow,GCSDur,ResTime, CalcEcoScoreDry, CalcEcoScoreWet, prevScore
  REAL:: TooWetStart
  REAL:: x0
  REAL,PARAMETER:: a=0.0, b=-1.0, lamda=1.0
  REAL, ALLOCATABLE :: TooWetTimeUp(:), TooWetCondUp(:), DryTimeUp(:), DryCondUp(:), DryTimeAdjUp(:), DryTimeLow(:), DryCondLow(:), DryTimeAdjLow(:)
  REAL, ALLOCATABLE :: WetTimeLow(:), WetCondLow(:), WetTimeAdjLow(:), TooWetTimeLow(:), TooWetCondLow(:), startScoreupExp(:), startScorelowExp(:)
  CHARACTER(20):: fmt
  CHARACTER(50)::infilename(12)
  
    
  TYPE InundType
    REAL :: CondScore, InundTime, InundTimeAdj
  END TYPE InundType
  TYPE (InundType), ALLOCATABLE::InundScore(:)

  TYPE ExpertERMDryType
    CHARACTER(50) :: ERMfileup, ERMfilelow
    INTEGER :: numERMup, numERMlow,numlinesup,numlineslow
    REAL :: TimeUp(12), CondUp(12), TimeLow(12), CondLow(12), TimeAdjUp(12), TimeAdjLow(12)
  END TYPE ExpertERMDryType
  TYPE (ExpertERMDryType), ALLOCATABLE :: ExpERMDry(:)

  TYPE ExpertERMWetType
    CHARACTER(50) :: ERMfileup, ERMfilelow
    INTEGER :: numERMup, numERMlow,numlinesup,numlineslow
    REAL :: TimeUp(12), CondUp(12), TimeLow(12), CondLow(12), TimeAdjUp(12), TimeAdjLow(12)
  END TYPE ExpertERMWetType
  TYPE (ExpertERMWetType), ALLOCATABLE :: ExpERMWet(:)

  TYPE ExpertERMTooWetGoodType
    CHARACTER(50) :: ERMfileup, ERMfilelow
    INTEGER :: numERMup, numERMlow,numlinesup,numlineslow
    REAL :: TimeUp(12), CondUp(12), TimeLow(12), CondLow(12), TimeAdjUp(12), TimeAdjLow(12)
  END TYPE ExpertERMTooWetGoodType
  TYPE (ExpertERMTooWetGoodType), ALLOCATABLE :: ExpERMTooWetG(:)

  TYPE ExpertERMTooWetPoorType
    CHARACTER(50) :: ERMfileup, ERMfilelow
    INTEGER :: numERMup, numERMlow,numlinesup,numlineslow
    REAL :: TimeUp(12), CondUp(12), TimeLow(12), CondLow(12), TimeAdjUp(12), TimeAdjLow(12)
  END TYPE ExpertERMTooWetPoorType
  TYPE (ExpertERMTooWetPoorType), ALLOCATABLE :: ExpERMTooWetP(:)

  
  CHARACTER(256) :: wd

  ALLOCATE (InundScore(8))
  ALLOCATE (TooWetTimeUp(4), TooWetCondUp(4), DryTimeUp(8), DryCondUp(8), DryTimeAdjUp(8), DryTimeLow(7), DryCondLow(7), DryTimeAdjLow(7))
  ALLOCATE (WetTimeLow(5), WetCondLow(5), WetTimeAdjLow(5), TooWetTimeLow(5), TooWetCondLow(5))
 
 !Open parameter files
 OPEN (unit=810,file='ConMod\ExpertERMs.txt', status='old')  !Contains file names containing expert response curves
 
 !Parameter files if all expert curves are combined into a single set of upper and lower response curves
 OPEN (unit=600,file='ConMod\WetUpperBound.txt',status='old')
 OPEN (unit=700,file='ConMod\TooWetUpperBound.txt', status='old')
 OPEN (unit=800,file='ConMod\DryUpperBound.txt', status='old')
 OPEN (unit=900,file='ConMod\DryLowerBound.txt', status='old')
 OPEN (unit=610,file='ConMod\WetLowerBound.txt', status='old')
 OPEN (unit=710,file='ConMod\TooWetLowerBound.txt', status='old')

 !Open output files
 OPEN (unit=200,file='ConMod\EcoScoreUp.txt',status='unknown')
 OPEN (unit=300,file='ConMod\EcoScoreLow.txt',status='unknown') 
 OPEN (unit=210,file='ConMod\EcoScoreUpExp.txt', status='unknown')
 OPEN (unit=310,file='ConMod\EcoScoreLowExp.txt', status='unknown')

 !----------------------------------------------------------------------
 !Read in expert ecological response curve parameters
 READ(810,*)numexperts
 READ(810,*)numexpertupper
 READ(810,*)numexpertlower
 READ(810,*)numexperttoowet

 ALLOCATE(ExpERMDry(numexperts))
 ALLOCATE(ExpERMWet(numexperts))
 ALLOCATE(ExpERMTooWetG(numexperttoowet))
 ALLOCATE(ExpERMTooWetP(numexperttoowet))
 ALLOCATE(startScoreupExp(numexpertupper))
 ALLOCATE(startScorelowExp(numexpertlower))


 DO i=1,numexpertupper
   READ(810,*)ExpERMDry(i)%numERMup
   READ(810,*)ExpERMDry(i)%ERMfileup
   OPEN (unit=i+50, file=ExpERMDry(i)%ERMfileup,status='old')
   READ (i+50,*)ExpERMDry(i)%numlinesup
   DO k=1,ExpERMDry(i)%numlinesup
     READ(i+50,*) ExpERMDry(i)%TimeUp(k), ExpERMDry(i)%CondUp(k)
     ExpERMDry(i)%TimeUp(k)=ExpERMDry(i)%TimeUp(k)*365 !Convert from years to days
   END DO

   READ(810,*)ExpERMWet(i)%ERMfileup
   OPEN (unit=i+50, file=ExpERMWet(i)%ERMfileup,status='old')
   READ (i+50,*)ExpERMWet(i)%numlinesup
   DO k=1,ExpERMWet(i)%numlinesup
     READ(i+50,*) ExpERMWet(i)%TimeUp(k), ExpERMWet(i)%CondUp(k)
     ExpERMWet(i)%TimeUp(k)=ExpERMWet(i)%TimeUp(k)*30.5 !Convert from months to days
   END DO
 END DO

 !Read in lower bounds for dry and wet conditions for each expert
 DO i=1,numexpertlower
   READ(810,*)ExpERMDry(i)%numERMlow
   READ(810,*)ExpERMDry(i)%ERMfilelow
   OPEN (unit=i+50, file=ExpERMDry(i)%ERMfilelow,status='old')
   READ (i+50,*)ExpERMDry(i)%numlineslow
   DO k=1,ExpERMDry(i)%numlineslow
     READ(i+50,*) ExpERMDry(i)%TimeLow(k), ExpERMDry(i)%CondLow(k)
     ExpERMDry(i)%TimeLow(k)=ExpERMDry(i)%TimeLow(k)*365 !Convert from years to days
   END DO

   READ(810,*)ExpERMWet(i)%ERMfilelow
   OPEN (unit=i+50, file=ExpERMWet(i)%ERMfilelow,status='old')
   READ (i+50,*)ExpERMWet(i)%numlineslow
   DO k=1,ExpERMWet(i)%numlineslow
     READ(i+50,*) ExpERMWet(i)%TimeLow(k), ExpERMWet(i)%CondLow(k)
     ExpERMWet(i)%TimeLow(k)=ExpERMWet(i)%TimeLow(k)*30.5 !Convert from years to days
   END DO
 END DO

 !Read in upper bounds for dry and wet conditions for each expert
 DO i=1,numexperttoowet  !Upper bound
   READ(810,*)ExpERMTooWetG(i)%ERMfileup
   OPEN(unit=i+60,file=ExpERMTooWetG(i)%ERMfileup,status='old')
   READ(i+60,*)ExpERMTooWetG(i)%numlinesup
   DO k=1,ExpERMTooWetG(i)%numlinesup
     READ(i+60,*)ExpERMTooWetG(i)%TimeUp(k), ExpERMTooWetG(i)%CondUp(k)
     ExpERMTooWetG(i)%TimeUp(k)=ExpERMTooWetG(i)%TimeUp(k)*30.5 !Convert from years to days
   END DO

   READ(810,*)ExpERMTooWetP(i)%ERMfileup
   OPEN(unit=i+60,file=ExpERMTooWetP(i)%ERMfileup,status='old')
   READ(i+60,*)ExpERMTooWetP(i)%numlinesup
   DO k=1,ExpERMTooWetP(i)%numlinesup
     READ(i+60,*)ExpERMTooWetP(i)%TimeUp(k), ExpERMTooWetP(i)%CondUp(k)
     ExpERMTooWetP(i)%TimeUp(k)=ExpERMTooWetP(i)%TimeUp(k)*30.5 !Convert from years to days
   END DO
 END DO

 !Read in lower bounds for too wet conditions (starting in good and poor condition) for each expert
 DO i=1,numexperttoowet  !Lower bound
   READ(810,*)ExpERMTooWetG(i)%ERMfilelow
   OPEN(unit=i+60,file=ExpERMTooWetG(i)%ERMfilelow,status='old')
   READ(i+60,*)ExpERMTooWetG(i)%numlineslow
   DO k=1,ExpERMTooWetG(i)%numlineslow
     READ(i+60,*)ExpERMTooWetG(i)%TimeLow(k), ExpERMTooWetG(i)%CondLow(k)
     ExpERMTooWetG(i)%TimeLow(k)=ExpERMTooWetG(i)%TimeLow(k)*30.5 !Convert from years to days
   END DO

   READ(810,*)ExpERMTooWetP(i)%ERMfilelow
   OPEN(unit=i+60,file=ExpERMTooWetP(i)%ERMfilelow,status='old')
   READ(i+60,*)ExpERMTooWetP(i)%numlineslow
   DO k=1,ExpERMTooWetP(i)%numlineslow
     READ(i+60,*)ExpERMTooWetP(i)%TimeLow(k), ExpERMTooWetP(i)%CondLow(k)
     ExpERMTooWetP(i)%TimeLow(k)=ExpERMTooWetP(i)%TimeLow(k)*30.5 !Convert from years to days
   END DO
 END DO

 !The following code reads in response curves for a pre-computed combined upper and lower bound for all experts.
 !Not currently used, but kept in case useful later. THis hasn't been recently checked 13/4/2014
 DO i=1, 8
   READ(600,*) InundScore(i)%InundTime, InundScore(i)%CondScore
   InundScore(i)%InundTime=InundScore(i)%InundTime*30.5  !Convert from months to days
   READ(800,*) DryTimeUp(i), DryCondUp(i)
   DryTimeUp(i)=DryTimeUp(i)*365  !Convert from years to days
 END DO

 DO i=1,7
   READ(900,*) DryTimeLow(i), DryCondLow(i)
   DryTimeLow(i)=DryTimeLow(i)*365  !Convert from years to days
 END DO

 DO i=1,4
   READ(700,*) TooWetTimeUp(i), TooWetCondUp(i)
   TooWetTimeUp(i)=TooWetTimeUp(i)*30.5
 END DO 

 DO i=1,5
   READ(610,*) WetTimeLow(i), WetCondLow(i)
   WetTimeLow(i)=WetTimeLow(i)*30.5  !Convert from months to days
 END DO

 DO i=1,5
   READ(710,*) TooWetTimeLow(i), TooWetCondLow(i)
   TooWetTimeLow(i)=TooWetTimeLow(i)*30.5
 END DO 

 !-----------------------------------------------------------------------------
  
  !Initialise variables  
  sumEcoScore=0.0
  sumEcoScoreup=0.0
  sumEcoScorelow=0.0
  startEcoScore=0.7
  ResTime=0

  
 test = GETCWD (wd)
 GWadj=0
 GWCount=0
 GWav=0
 GWEndCount=0
 GWCountAbove9=0
 TooWetCountUp(6)=0
 TooWetCountLow(6)=0

 CALL CompareInundationArea (x0, maxk, GWCode, GWDepthCode)  !Calculate inundation duration of the GCS

!--------------------------------------------------------------------------------
!ECOLOGICAL RESPONSE MODEL



!--------------------------------------------------------------------------------


 i=1

   DO k=1,maxk   !loop through each event

      !-------------------------------------------------------------------
      !DRY PERIOD
      !------------------------------------------------------------------
      IF (GCSevent(k)%classif==1) THEN !dry period
       GWCount=0
       GWadj=0.0
       GWav=0.0
       GWEndCount=0
       GWCountAbove9=0
       DryGW=0

       DO e=1,numexpertupper
         DryGWExp(e)=0
         WetGWExp(e)=0
       END DO

       DO e=1,numexpertlower
         DryGWlowExp(e)=0
       END DO

       DO m=1,GCSevent(k)%Dur    !loop through each day in the event

        !Calculate GW % scale factor (Specify GWDepth 1 or 2)
        IF (GWDepthCode==1) THEN
          obj(i)%GWDepth=obj(i)%GWDepth   !Use observed GW levels
        ELSE
          obj(i)%GWDepth=obj(i)%GWDepth2  !Use halved observed GW levels
        END IF

        GWscale=1/(1+exp(-lamda*(obj(i)%GWDepth-x0)))   !Estimate the percentage access to groundwater based on RRG root depths


         IF (k==1 .and. m==1) THEN                        !Set initial condition score
           startScoreup=startEcoScore
           startScorelow=startEcoScore
           DO e=1,numexperts
             startScoreupExp(e)=startEcoScore
             startScorelowExp(e)=startEcoScore
           END DO
         ELSE IF (m==1) THEN
           startScoreup= obj(i-1)%ecoScoreup  
           startScorelow= obj(i-1)%ecoScorelow
           DO e=1,numexperts
             startScoreupExp(e)=obj(i-1)%ecoScoreupExp(e)
             startScorelowExp(e)=obj(i-1)%ecoScorelowExp(e)
           END DO
         END IF  
         
         IF (m==1) THEN
           !Calculate equivalent duration based on starting condition for upper bound (as lookup table assumes score=1.0)
         !   startScoreup=0.0
           DurAdj=DurAdjustDry(DryTimeUp(:),DryCondUp(:),startScoreup,8)

           ! Calculate adjusted time in lookup take to account for a different start score
           DO t=1,8
             DryTimeAdjUp(t)=DryTimeUp(t)-DurAdj
           END DO

           !Calculate adjusted time for individual experts
           DO e=1, numexpertupper
             DurAdj=DurAdjustDry(ExpERMDry(e)%TimeUp(:), ExpERMDry(e)%CondUp(:),startScoreupExp(e),ExpERMDry(e)%numlinesup)
             DO t=1,ExpERMDry(e)%numlinesup
               ExpERMDry(e)%TimeAdjUp(t)=ExpERMDry(e)%TimeUp(t)-DurAdj
             END DO
           END DO  !End from e=1 to numexpertupper

           !Calculate equivalent duration based on starting condition for lower bound (as lookup table assumes score=1.0)
           DurAdj=DurAdjustDry(DryTimeLow(:),DryCondLow(:),startScorelow,7)          

           ! Calculate adjusted time in lookup take to account for a different start score
           DO t=1,7
             DryTimeAdjLow(t)=DryTimeLow(t)-DurAdj
           END DO

           !Calculate adjusted time for individual experts
           DO e=1, numexpertlower
             DurAdj=DurAdjustDry(ExpERMDry(e)%TimeLow(:), ExpERMDry(e)%CondLow(:),startScorelowExp(e),ExpERMDry(e)%numlineslow)
             DO t=1,ExpERMDry(e)%numlineslow
               ExpERMDry(e)%TimeAdjLow(t)=ExpERMDry(e)%TimeLow(t)-DurAdj
             END DO
           END DO  !End from e=1 to numexpertlower
           
         END IF !end if m==1 loop

         IF (GWscale>0.01) THEN
           GWCount=GWCount+1                  !Count the number of time steps where % access to GW is >0.01
           GWadj=GWadj+GWscale
           GWav=GWadj/REAL(GWCount)           !Calculate the average % access over this duration to use in adjusting condition scores
         END IF

        !-------------------------------------------------- 
            !Calculate current score for upper bound.
        !----------------------------------------------

           !Determine if was in previously too wet condition, in which case an initial improvement may be observed
           DO e=1,numexpertupper
             IF (TooWetCountUp(e)>0) THEN  ! was previously too wet
               IF (e==1) THEN  !Expert 1
                 IF (m<=183)  THEN
                   obj(i)%ecoScoreupExp(e)=min(obj(i-1)%ecoScoreupExp(e) + (0.1/183.0),1.0) ! score improves for first 6 months at a rate of 0.1 per 6 months to a maximum of 1
                 ELSE IF (m<365) THEN
                   obj(i)%ecoScoreupExp(e)=obj(i-1)%ecoScoreupExp(e)  !Score stabilises until 1 year has passed
                 ELSE
                   TooWetCountUp(e)=0  !Reset to zero so will follow normal drought trajectory
                 END IF
               END IF
               IF (e==2) THEN   !Expert 2
                 IF (startScoreupExp(e)>0.4) THEN
                   IF (m<1095) THEN  !Condition stays the same for an upper bound of 3 years
                     obj(i)%ecoScoreupExp(e)=obj(i-1)%ecoScoreupExp(e)
                   ELSE
                     TooWetCountUp(e)=0
                   END IF
                 ELSE     !If score<=0.4
                   IF (m<1095) THEN  !Gradual improvement over an upper bound of 3 years - amount estimated
                     obj(i)%ecoScoreupExp(e)=min(obj(i-1)%ecoScoreupExp(e) + (0.2/1095),1.0)
                   ELSE
                     TooWetCountUp(e)=0
                   END IF
                 END IF
               END IF
               IF (e==3 .or. e==4) THEN   !Expert 3 and 4
                 IF (m<=183) THEN
                   obj(i)%ecoScoreupExp(e)=min(obj(i-1)%ecoScoreupExp(e) + (0.1/183.0),1.0)  !Score improves 0.1 over 6 months
                 ELSE IF (m<=730) THEN
                   obj(i)%ecoScoreupExp(e)=obj(i-1)%ecoScoreupExp(e)   !Condition stabilises up to 2 years
                 ELSE
                   TooWetCountUp(e)=0
                 END IF
               END IF
             END IF !End if toowetcount>0
           END DO !End from 1 to numexpertupper
           
           
           IF (GWCode==1 .and. GWscale>0.9) THEN   !Assume follows curve for wet period upper bound
             GWCountAbove9=GWCountAbove9+1
             IF (GWCountAbove9==1) THEN
             !Calculate duration adjustment
              
               IF (i==1) THEN
                 prevScore=startScoreup
               ELSE
                 prevScore=obj(i-1)%ecoScoreup
               END IF
               
               DurAdj=DurAdjustWet(InundScore(:)%InundTime,InundScore(:)%CondScore,prevScore,8)

               ! Calculate adjusted time in lookup take to account for a different start score
               DO t=1,8
                 InundScore(t)%InundTimeAdj=InundScore(t)%InundTime-DurAdj
               END DO
             END IF   !End if GWCountAbove9==1

               !Calculate for all experts
                
             DO e=1, numexpertupper
               IF (TooWetCountUp(e)==0) THEN             !Only consider GW if isn't following a too wet period, as condition has already been calculated above
                 WetGWExp(e)=WetGWExp(e)+1               !Count number of days GW above 90% access for each expert
                 IF (WetGWExp(e)==1) THEN
                   IF(i==1) THEN
                     prevScore=startScoreupExp(e)
                   ELSE
                     prevScore=obj(i-1)%ecoScoreupExp(e)
                   END IF
                   DurAdj=DurAdjustWet(ExpERMWet(e)%TimeUp(:), ExpERMWet(e)%CondUp(:),prevScore,ExpERMWet(e)%numlinesup)
                   DO t=1,ExpERMWet(e)%numlinesup
                     ExpERMWet(e)%TimeAdjUp(t)=ExpERMWet(e)%TimeUp(t)-DurAdj
                   END DO
                 END IF  !If WetGWExp(e)==1
               END IF   ! If TooWetCountUp(e)==0
             END DO  !End from e=1 to numexpertupper
               
             
             !Calculate score for if GWscale>0.9
             IF (i==1) THEN
               prevScore=startScoreup
             ELSE
               prevScore=obj(i-1)%ecoScoreup
             END IF
             obj(i)%ecoScoreup=CalcEcoScoreWet(InundScore(:)%InundTimeAdj,InundScore(:)%CondScore,prevScore,8,GWCountAbove9)

             !Calculate for all experts
             DO e=1,numexpertupper
             IF (TooWetCountUp(e)==0) THEN
               IF(i==1) THEN
                 prevScore=startScoreupExp(e)
               ELSE
                 prevScore=obj(i-1)%ecoScoreupExp(e)
               END IF
               obj(i)%ecoScoreupExp(e)=CalcEcoScoreWet(ExpERMWet(e)%TimeAdjUp(:),ExpERMWet(e)%CondUp(:),prevScore,ExpERMWet(e)%numlinesup,WetGWExp(e))
             END IF
             END DO
               
           ELSE IF (GWCode==1 .and. GWscale>0.01) THEN   !Else if GW<=0.9)
             GWCountAbove9=0
             DryGW=DryGW+1   !Start dry period duration 
             IF (DryGW==1) THEN    !Adjust durations based on current score
               IF(i==1) THEN
                 prevScore=startScoreup
               ELSE
                 prevScore=obj(i-1)%ecoScoreup
               END IF
               DurAdj=DurAdjustDry(DryTimeUp(:),DryCondUp(:),prevScore,8)

               ! Calculate adjusted time in lookup take to account for a different start score
               DO t=1,8
                 DryTimeAdjUp(t)=DryTimeUp(t)-DurAdj
               END DO
             END IF

               !Calculate adjusted time for individual experts
             DO e=1, numexpertupper
               IF (TooWetCountUp(e)==0) THEN              !Only calculate if not following a wet period, as this has already been calculated (or if dry period has been sufficiently long)
                 DryGWExp(e)=DryGWExp(e)+1
                   IF (DryGWExp(e)==1) THEN
                     IF (i==1) THEN
                       prevScore=startScoreupExp(e)
                     ELSE
                       prevScore=obj(i-1)%ecoScoreupExp(e)
                     END IF               
                    
                     DurAdj=DurAdjustDry(ExpERMDry(e)%TimeUp(:), ExpERMDry(e)%CondUp(:),prevScore,ExpERMDry(e)%numlinesup)
   
                     DO t=1,ExpERMDry(e)%numlinesup 
                       ExpERMDry(e)%TimeAdjUp(t)=ExpERMDry(e)%TimeUp(t)-DurAdj
                     END DO
                   END IF
                END IF
             END DO  !End from e=1 to numexpertupper
           

             DO t=1,7
               DryTimeAdjUp(t)=DryTimeAdjUp(t)+GWscale
             END DO
            
             IF (i==1) THEN
               prevScore=startScoreup
             ELSE
               prevScore=obj(i-1)%ecoScoreup
             END IF
             obj(i)%ecoScoreup=CalcEcoScoreDry(DryTimeAdjUp(:),DryCondUp(:),prevScore,8,DryGW)

             !Calculate for all experts
             DO e=1,numexpertupper
             IF (TooWetCountUp(e)==0) THEN
               DO t=1,ExpERMDry(e)%numlinesup
                 ExpERMDry(e)%TimeAdjUp(t)=ExpERMDry(e)%TimeAdjUp(t)+GWscale  !By adding the % GW access (as a factor), it takes slightly longer for the condition to decline than without GW access
               END DO
               IF(i==1) THEN
                 prevScore=startScoreupExp(e)
               ELSE
                 prevScore=obj(i-1)%ecoScoreupExp(e)
               END IF
              obj(i)%ecoScoreupExp(e)=CalcEcoScoreDry(ExpERMDry(e)%TimeAdjUp(:),ExpERMDry(e)%CondUp(:),prevScore,ExpERMDry(e)%numlinesup,DryGWExp(e))
             END IF
             END DO

           ELSE   !If GW<0.01               !Follows normal dry period curve
             GWCountAbove9=0
             DryGW=DryGW+1           
             IF (i==1) THEN
               prevScore=startScoreup
             ELSE
               prevScore=obj(i-1)%ecoScoreup
             END IF

             IF (DryGW==1) THEN    !Adjust durations based on current score
               IF(i==1) THEN
                 prevScore=startScoreup
               ELSE
                 prevScore=obj(i-1)%ecoScoreup
               END IF
               DurAdj=DurAdjustDry(DryTimeUp(:),DryCondUp(:),prevScore,8)

               ! Calculate adjusted time in lookup take to account for a different start score
               DO t=1,8
                 DryTimeAdjUp(t)=DryTimeUp(t)-DurAdj
               END DO
             END IF

               !Calculate adjusted time for individual experts
               DO e=1, numexpertupper
                 IF (TooWetCountUp(e)==0) THEN
                   DryGWExp(e)=DryGWExp(e)+1
                   IF (DryGWExp(e)==1) THEN
                     IF (i==1) THEN
                       prevScore=startScoreupExp(e)
                     ELSE
                       prevScore=obj(i-1)%ecoScoreupExp(e)
                     END IF
                     DurAdj=DurAdjustDry(ExpERMDry(e)%TimeUp(:), ExpERMDry(e)%CondUp(:),startScoreupExp(e),ExpERMDry(e)%numlinesup)
                     DO t=1,ExpERMDry(e)%numlinesup
                       ExpERMDry(e)%TimeAdjUp(t)=ExpERMDry(e)%TimeUp(t)-DurAdj
                     END DO
                   END IF
                 END IF
               END DO  !End from e=1 to numexpertupper
            
             obj(i)%ecoScoreup=CalcEcoScoreDry(DryTimeAdjUp(:),DryCondUp(:),prevScore,8,DryGW)

             !Calculate for all experts
             DO e=1,numexpertupper
             IF (TooWetCountUp(e)==0) THEN
               IF(i==1) THEN
                 prevScore=startScoreupExp(e)
               ELSE
                 prevScore=obj(i-1)%ecoScoreupExp(e)
               END IF
               obj(i)%ecoScoreupExp(e)=CalcEcoScoreDry(ExpERMDry(e)%TimeAdjUp(:),ExpERMDry(e)%CondUp(:),prevScore,ExpERMDry(e)%numlinesup,DryGWExp(e))
             END IF
             END DO
           
           END IF   !End if GW>0.9 or <0.01

            !-------------------------------------------------------
            !Calculate current score for LOWER bound. DRY CONDITIONS
            !-------------------------------------------------------

            !Determine if was in previously too wet condition, in which case an initial improvement may be observed
           DO e=1,numexpertlower
             IF (TooWetCountLow(e)>0) THEN  ! was previously too wet
               IF (e==1) THEN  
                 IF (m<=183)  THEN
                   obj(i)%ecoScorelowExp(e)=min(1.0,obj(i-1)%ecoScorelowExp(e) + (0.1/183.0)) ! score improves for first 6 months at a rate of 0.1 per 6 months
                 ELSE IF (m<365) THEN
                   obj(i)%ecoScorelowExp(e)=obj(i-1)%ecoScorelowExp(e)  !Score stabilises until 1 year has passed
                 ELSE
                   TooWetCountLow(e)=0  !Reset to zero so will follow normal drought trajectory
                 END IF
               END IF
               IF (e==2) THEN   
                 IF (startScorelowExp(e)>0.4) THEN
                   IF (m<730) THEN  !Condition stays the same for an lower bound of 2 years
                     obj(i)%ecoScorelowExp(e)=obj(i-1)%ecoScorelowExp(e)
                   ELSE
                     TooWetCountLow(e)=0
                   END IF
                 ELSE     !If score<=0.4
                   IF (m<730) THEN  !Gradual improvement over an upper bound of 3 years - amount estimated
                     obj(i)%ecoScorelowExp(e)=min(1.0,obj(i-1)%ecoScorelowExp(e) + (0.1/730))
                   ELSE
                     TooWetCountLow(e)=0
                   END IF
                 END IF
               END IF
               IF (e==3) THEN   
                 IF (m<=183) THEN
                   obj(i)%ecoScorelowExp(e)=min(1.0,obj(i-1)%ecoScorelowExp(e) + (0.1/183.0))  !Score improves 0.1 over 6 months
                 ELSE IF (m<=365) THEN
                   obj(i)%ecoScorelowExp(e)=obj(i-1)%ecoScorelowExp(e)   !Condition stabilises up to 1 year
                 ELSE
                   TooWetCountLow(e)=0
                 END IF
               END IF
               IF (e==4) THEN  
                 IF (startScorelowExp(e)>0.4) THEN
                   IF (m<=183) THEN
                     obj(i)%ecoScorelowExp(e)=min(1.0, obj(i-1)%ecoScorelowExp(e) + (0.1/183.0))  !Same as upper bound
                   ELSE IF (m<=730) THEN
                     obj(i)%ecoScorelowExp(e)=obj(i-1)%ecoScorelowExp(e)
                   ELSE
                     TooWetCountLow(e)=0
                   END IF
                 ELSE  !start score <=0.4
                   obj(i)%ecoScorelowExp(e)=max(0.0,obj(i-1)%ecoScorelowExp(e) - (startScorelowExp(e)/730)) !Reduces to 0 after 2 years
                 END IF
               END IF
             END IF !End if toowetcount>0
           END DO !End from 1 to numexpertlower


               DO e=1,numexpertlower
                 IF (TooWetCountLow(e)==0) THEN   !Determine if start of dry period, or still recovering from too wet
                   DryGWlowExp(e)=DryGWlowExp(e)+1
                   IF (DryGWlowExp(e)==1) THEN    !start of dry period, calculate modified duration
                     IF (i==1) THEN 
                       prevScore=startScorelowExp(e)
                     ELSE
                       prevScore=obj(i-1)%ecoScorelowExp(e)
                     END IF
                     DurAdj=DurAdjustDry(ExpERMDry(e)%TimeLow(:), ExpERMDry(e)%CondLow(:),prevScore,ExpERMDry(e)%numlineslow)
                     DO t=1,ExpERMDry(e)%numlineslow
                       ExpERMDry(e)%TimeAdjLow(t)=ExpERMDry(e)%TimeLow(t)-DurAdj
                     END DO
                   END IF
                 END IF
               END DO


               IF (GWCode==1) THEN
                 IF(GWscale>0.01) THEN
                   DO t=1,7
                     DryTimeAdjLow(t)=DryTimeAdjLow(t)+GWscale
                   END DO
                   DO e=1,numexpertlower
                     IF (TooWetCountLow(e)==0) THEN
                       DO t=1,ExpERMDry(e)%numlinesup
                         ExpERMDry(e)%TimeAdjLow(t)=ExpERMDry(e)%TimeAdjLow(t)+GWscale
                       END DO
                     END IF
                   END DO
                 ELSE
                   IF (GWCount>1) THEN
                     GWEndCount=GWEndCount+1
                     IF (GWEndCount<=GWCount) THEN
                         DO t=1,7
                           DryTimeAdjLow(t)=DryTimeAdjLow(t)-GWav
                         END DO
                         DO e=1,numexpertlower
                           DO t=1,ExpERMDry(e)%numlinesup
                             ExpERMDry(e)%TimeAdjLow(t)=ExpERMDry(e)%TimeAdjLow(t)-GWav
                           END DO
                         END DO                                                           
                     END IF   !If GWEndCount<=GWCount
                   END IF     !If GWCount>1
                 END IF    !If GWscale>0.01
               END IF    !If GWCode==1

           IF (i==1) THEN
             prevScore=startScorelow
           ELSE
             prevScore=obj(i-1)%ecoScorelow
           END IF
           
           obj(i)%ecoScorelow=CalcEcoScoreDry(DryTimeAdjLow(:),DryCondLow(:),prevScore,7,m)
           
           !For all experts
           DO e=1,numexpertlower
             IF (TooWetCountLow(e)==0) THEN
               IF (i==1) THEN
                 prevScore=startScorelowExp(e)
               ELSE
                 prevScore=obj(i-1)%ecoScorelowExp(e)
               END IF
               obj(i)%ecoScorelowExp(e)=CalcEcoScoreDry(ExpERMDry(e)%TimeAdjLow(:),ExpERMDry(e)%CondLow(:),prevScore,ExpERMDry(e)%numlineslow,DryGWlowExp(e))
             END IF
           END DO

         IF(i>1) THEN
           IF (obj(i-1)%ecoScoreup<0.00001) THEN
             obj(i)%ecoScoreup=0.0
           END IF
           DO e=1,numexpertupper 
             IF (obj(i-1)%ecoScoreupExp(e)<0.00001) THEN
               obj(i)%ecoScoreupExp(e)=0.0
             END IF
           END DO 
           IF (obj(i-1)%ecoScorelow<0.00001) THEN
             obj(i)%ecoScorelow=0.0
           END IF 
           DO e=1,numexpertlower 
             IF (obj(i-1)%ecoScorelowExp(e)<0.00001) THEN
               obj(i)%ecoScorelowExp(e)=0.0
             END IF
           END DO 
         END IF
         

  !       WRITE(200,*)obj(i)%ecoScoreup            !Only use if calculating aggregated upper and lower bounds for all experts
  !       WRITE(300,*)obj(i)%ecoScorelow
         fmt = '(5(3x,f8.6))'
         WRITE(210,fmt)(obj(i)%ecoScoreupExp(e),e=1,numexpertupper)
         WRITE(310,fmt)(obj(i)%ecoScorelowExp(e),e=1,numexpertlower)
                          
         i=i+1 
        END DO
     ELSE 
          !------------------------------------------------------------------------------------
          !WET PERIOD 


          !--------------------------------------------------------------------------------
    
       DO m=1,GCSevent(k)%Dur  
         IF (k==1 .and. m==1) THEN
           startScoreup=startEcoScore
           startScorelow=startEcoScore
           DO e=1,numexperts
             startScoreupExp(e)=startEcoScore
             startScorelowExp(e)=startEcoScore
           END DO
         ELSE IF (m==1) THEN  !Start of the event - calculate initial score and duration adjustment
           startScoreup= obj(i-1)%ecoScoreup  
           startScorelow= obj(i-1)%ecoScorelow
           DO e=1,numexperts
             startScoreupExp(e)=obj(i-1)%ecoScoreupExp(e)
             startScorelowExp(e)=obj(i-1)%ecoScorelowExp(e)
           END DO
         END IF
         
         IF (m==1) THEN
           !Calculate equivalent duration based on starting condition for upper bound (as lookup table assumes score=0.2)
           DurAdj=DurAdjustWet(InundScore(:)%InundTime,InundScore(:)%CondScore,startScoreup,8)

           ! Calculate adjusted time in lookup take to account for a different start score
           DO t=1,8
             InundScore(t)%InundTimeAdj=InundScore(t)%InundTime-DurAdj
           END DO

           !Calculate adjusted time for individual experts
           DO e=1, numexpertupper
             DurAdj=DurAdjustWet(ExpERMWet(e)%TimeUp(:), ExpERMWet(e)%CondUp(:),startScoreupExp(e),ExpERMWet(e)%numlinesup)
             DO t=1,ExpERMWet(e)%numlinesup
               IF (e==2) THEN   !If expert = 2
                 ExpERMWet(e)%TimeAdjUp(t)=ExpERMWet(e)%TimeUp(t)-DurAdj + 30   !Add a 1 month lag to the start of response
               ELSE
                 ExpERMWet(e)%TimeAdjUp(t)=ExpERMWet(e)%TimeUp(t)-DurAdj
               END IF
             END DO
           END DO  !End from e=1 to numexpertupper

        
           !Calculate equivalent duration based on starting condition for lower bound (as lookup table assumes score=0.2)
           DurAdj=DurAdjustWet(WetTimeLow(:),WetCondLow(:),startScorelow,5)

           ! Calculate adjusted time in lookup take to account for a different start score
           DO t=1,5
             WetTimeAdjLow(t)=WetTimeLow(t)-DurAdj + 91  !Add 91 to account for 3 month initial lag   
           END DO  

           !Calculate adjusted time for individual experts
           DO e=1, numexpertlower
             DurAdj=DurAdjustWet(ExpERMWet(e)%TimeLow(:), ExpERMWet(e)%CondLow(:),startScorelowExp(e),ExpERMWet(e)%numlineslow)
             DO t=1,ExpERMWet(e)%numlineslow
               IF (e==2) THEN  !If expert = 2
                 ExpERMWet(e)%TimeAdjLow(t)=ExpERMWet(e)%TimeLow(t)-DurAdj + 91   !Add 91 to account for 3 month initial lag 
               ELSE
                 ExpERMWet(e)%TimeAdjLow(t)=ExpERMWet(e)%TimeLow(t)-DurAdj
               END IF
             END DO
           END DO  !End from e=1 to numexpertlower
        
         END IF  !End if m==1 loop
    
            !-----------------------------------------------  
            !Calculate current score for upper bound
            !----------------------------------------------
         IF (m<TooWetTimeUp(1)) THEN    !Check if duration of event is less than threshold for too wet
           IF (i==1) THEN
             prevScore=startScoreup
           ELSE
             prevScore=obj(i-1)%ecoScoreup
           END IF
           
           obj(i)%ecoScoreup=CalcEcoScoreWet(InundScore(:)%InundTimeAdj,InundScore(:)%CondScore,prevScore,8,m)
         ELSE
           IF (m==732) THEN       !24 months
             TooWetCondUp(1)= obj(i-1)%ecoScoreup
             PLInup=TooWetTimeUp(2)
             PLConUp=TooWetCondUp(2)
             PLInLow=TooWetTimeUp(1)
             PLConLow=TooWetCondUp(1)
                 
             obj(i)%ecoScoreup=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
       !      EXIT 
           ELSE              
             DO t=1,4
               IF (TooWetTimeUp(t)>=m) THEN
                 PLInUp=TooWetTimeUp(t)
                 PLConUp=TooWetCondUp(t)
                 PLInLow=TooWetTimeUp(t-1)
                 PLConLow=TooWetCondUp(t-1)

                 obj(i)%ecoScoreup=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
                 EXIT
               END IF  !End if too wet time>m
             END DO
           END IF !End if m==732
         END IF ! End if m< too wet time up(1)
     
         !Calculate current upper score for each expert   
         DO e=1,numexpertupper
           IF(i==1) THEN
             prevScore=startScoreupExp(e)
           ELSE
             prevScore=obj(i-1)%ecoScoreupExp(e)
           END IF
           IF (e==5) THEN  !Landholder 2 - no too wet defined so just duration of 30,000 and follow normal wet curve
             TooWetStart=30000
           ELSE
             IF(startScoreupExp(e)>0.4) THEN
               TooWetStart=ExpERMTooWetG(e)%TimeUp(1)   !Good condition
             ELSE
               TooWetStart=ExpERMTooWetP(e)%TimeUp(1)   ! poor condition
             END IF
           END IF  !End if e==5
           IF (m<TooWetStart) THEN
             TooWetCountUp(e)=0
             IF (e==2) THEN  !If Expert 2
               IF (m<=30) THEN  !lag for 1 month
                 IF (i==1) THEN
                   obj(i)%ecoScoreupExp(e)=startScoreupExp(e)
                 ELSE
                   obj(i)%ecoScoreupExp(e)=obj(i-1)%ecoScoreupExp(e)
                 END IF  !End if i==1
               ELSE
                  obj(i)%ecoScoreupExp(e)=CalcEcoScoreWet(ExpERMWet(e)%TimeAdjUp(:),ExpERMWet(e)%CondUp(:),prevScore,ExpERMWet(e)%numlinesup,m)
               END IF  !End if m<=30
             ELSE  !If e not = 2
               obj(i)%ecoScoreupExp(e)=CalcEcoScoreWet(ExpERMWet(e)%TimeAdjUp(:),ExpERMWet(e)%CondUp(:),prevScore,ExpERMWet(e)%numlinesup,m)
             END IF  !End if e==2
           ELSE   !If m>TooWetStart        
             TooWetCountUp(e)=TooWetCountUp(e)+1
             IF(TooWetCountUp(e)==1) THEN  !Calculate current score at start of too wet period
               IF (startScoreupExp(e)>0.4) THEN
                 ExpERMTooWetG(e)%CondUp(1)=obj(i-1)%ecoScoreupExp(e)
                 PLInup=ExpERMTooWetG(e)%TimeUp(2)
                 PLConUp=ExpERMTooWetG(e)%CondUp(2)
                 PLInLow=ExpERMTooWetG(e)%TimeUp(1)
                 PLConLow=ExpERMTooWetG(e)%CondUp(1)
                 
                 obj(i)%ecoScoreupExp(e)=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow

                 IF (e==4) THEN  !Modify time of death for Expert 4 based on extrapolating previous points
                   ExpERMTooWetG(e)%TimeUp(3)=PLInup-(PLConUp-ExpERMTooWetG(e)%CondUp(3))*(PLInLow-PLInup)/(PLConLow-PLConUp)
                 END IF
          !       EXIT 
               ELSE  !If startScoreupExp(e) <=0.4
                 ExpERMTooWetP(e)%CondUp(1)=obj(i-1)%ecoScoreupExp(e)
                 PLInup=ExpERMTooWetP(e)%TimeUp(2)
                 PLConUp=ExpERMTooWetP(e)%CondUp(2)
                 PLInLow=ExpERMTooWetP(e)%TimeUp(1)
                 PLConLow=ExpERMTooWetP(e)%CondUp(1)
                 
                 obj(i)%ecoScoreupExp(e)=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow

                 IF (e==4) THEN  !Modify time of death for Expert 4 based on extrapolating previous points
                   ExpERMTooWetP(e)%TimeUp(3)=PLInup-(PLConUp-ExpERMTooWetP(e)%CondUp(3))*(PLInLow-PLInup)/(PLConLow-PLConUp)
                 END IF
             !    EXIT
               END IF  !End if prev score >0.4
             ELSE  !if TooWetCount>1
               IF (startScoreupExp(e)>0.4) THEN
         tloop:  DO t=1,ExpERMTooWetG(e)%numlinesup
                   IF (ExpERMTooWetG(e)%TimeUp(t)>=m) THEN
                     PLInUp=ExpERMTooWetG(e)%TimeUp(t)
                     PLConUp=ExpERMTooWetG(e)%CondUp(t)
                     PLInLow=ExpERMTooWetG(e)%TimeUp(t-1)
                     PLConLow=ExpERMTooWetG(e)%CondUp(t-1)

                     obj(i)%ecoScoreupExp(e)=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
                     EXIT tloop
                   END IF  !End if too wet time>m
                 END DO tloop
               ELSE !If startScoreupExp(e)<=0.4
          tloop2: DO t=1,ExpERMTooWetP(e)%numlinesup
                   IF (ExpERMTooWetP(e)%TimeUp(t)>=m) THEN
                     PLInUp=ExpERMTooWetP(e)%TimeUp(t)
                     PLConUp=ExpERMTooWetP(e)%CondUp(t)
                     PLInLow=ExpERMTooWetP(e)%TimeUp(t-1)
                     PLConLow=ExpERMTooWetP(e)%CondUp(t-1)

                     obj(i)%ecoScoreupExp(e)=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
                     EXIT tloop2
                   END IF  !End if too wet time>m
                 END DO tloop2
               END IF !End if prevScore>0.4
             END IF !End if TooWetCount>1
           END IF ! End if m< too wet start
         END DO

            !----------------------------------------------------
            !Calculate current score for lower bound
            !----------------------------------------------------
          
             IF (m<TooWetTimeLow(1)) THEN    !Check if duration of event is less than threshold for too wet
               IF (m<=91) THEN   !No change in condition for first 3 months
                 IF (i==1) THEN
                   obj(i)%ecoScorelow=startScorelow
                 ELSE
                   obj(i)%ecoScorelow=obj(i-1)%ecoScorelow
                 END IF
               ELSE   !if m>91
                 IF (i==1) THEN
                   prevScore=startScorelow
                 ELSE
                   prevScore=obj(i-1)%ecoScorelow
                 END IF
           
                 obj(i)%ecoScorelow=CalcEcoScoreWet(WetTimeAdjLow(:),WetCondLow(:),prevScore,5,m)
               END IF  !End if m<=91
             ELSE
               IF (m==214) THEN       !7 months
                 TooWetCondLow(1)= obj(i-1)%ecoScorelow 
               END IF 
                            
               DO t=1,5
                 IF (TooWetTimeLow(t)>=m) THEN
                   PLInUp=TooWetTimeLow(t)
                   PLConUp=TooWetCondLow(t)
                   PLInLow=TooWetTimeLow(t-1)
                   PLConLow=TooWetCondLow(t-1)

                   obj(i)%ecoScorelow=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
                   
                   EXIT
                 END IF
               END DO
             END IF     
          

             !Calculate current lower score for each expert
             DO e=1,numexpertlower
               IF(i==1) THEN
                 prevScore=startScorelowExp(e)
               ELSE
                 prevScore=obj(i-1)%ecoScorelowExp(e)
               END IF
               IF (e==5) THEN  !Landholder 2 - no too wet defined so just duration of 30,000 and follow normal wet curve
                 TooWetStart=30000
               ELSE
                 IF(startScorelowExp(e)>0.4) THEN
                   TooWetStart=ExpERMTooWetG(e)%TimeLow(1)   !Good condition
                 ELSE
                   TooWetStart=ExpERMTooWetP(e)%TimeLow(1)   ! poor condition
                 END IF
               END IF  !End if e==5
               
               IF (m<TooWetStart) THEN
                 TooWetCountLow(e)=0
                 IF (e==2) THEN  !If Expert 2
                   IF (m<=91) THEN  !lag for 3 months
                     IF (i==1) THEN
                       obj(i)%ecoScorelowExp(e)=startScorelowExp(e)
                     ELSE
                       obj(i)%ecoScorelowExp(e)=obj(i-1)%ecoScorelowExp(e)
                     END IF  !End if i==1
                   ELSE
                     obj(i)%ecoScorelowExp(e)=CalcEcoScoreWet(ExpERMWet(e)%TimeAdjLow(:),ExpERMWet(e)%CondLow(:),prevScore,ExpERMWet(e)%numlineslow,m)
                   END IF  !End if m<=91
                 ELSE  !If e not = 2
                   obj(i)%ecoScorelowExp(e)=CalcEcoScoreWet(ExpERMWet(e)%TimeAdjLow(:),ExpERMWet(e)%CondLow(:),prevScore,ExpERMWet(e)%numlineslow,m)
                 END IF  !End if e==2
               ELSE   !If m>TooWetStart
                 TooWetCountLow(e)=TooWetCountLow(e)+1
                 IF(TooWetCountLow(e)==1) THEN  !Calculate current score at start of too wet period
                   IF (startScorelowExp(e)>0.4) THEN
                     ExpERMTooWetG(e)%CondLow(1)=obj(i-1)%ecoScorelowExp(e)
                     PLInup=ExpERMTooWetG(e)%TimeLow(2)
                     PLConUp=ExpERMTooWetG(e)%CondLow(2)
                     PLInLow=ExpERMTooWetG(e)%TimeLow(1)
                     PLConLow=ExpERMTooWetG(e)%CondLow(1)
                 
                     obj(i)%ecoScorelowExp(e)=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow

                     IF (e==4) THEN  !Modify time of death for Expert 4 based on extrapolating previous points
                       ExpERMTooWetG(e)%TimeLow(3)=PLInup-(PLConUp-ExpERMTooWetG(e)%CondLow(3))*(PLInLow-PLInup)/(PLConLow-PLConUp)
                     END IF
                   ELSE  !If startScorelowExp(e) <=0.4
                     ExpERMTooWetP(e)%CondLow(1)=obj(i-1)%ecoScorelowExp(e)
                     PLInup=ExpERMTooWetP(e)%TimeLow(2)
                     PLConUp=ExpERMTooWetP(e)%CondLow(2)
                     PLInLow=ExpERMTooWetP(e)%TimeLow(1)
                     PLConLow=ExpERMTooWetP(e)%CondLow(1)
                 
                     obj(i)%ecoScorelowExp(e)=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
                   END IF  !End if prev score >0.4
                 ELSE  !if TooWetCount>1
                   IF (startScorelowExp(e)>0.4) THEN
            tloop3:  DO t=1,ExpERMTooWetG(e)%numlineslow
                       IF (ExpERMTooWetG(e)%TimeLow(t)>=m) THEN
                         PLInUp=ExpERMTooWetG(e)%TimeLow(t)
                         PLConUp=ExpERMTooWetG(e)%CondLow(t)
                         PLInLow=ExpERMTooWetG(e)%TimeLow(t-1)
                         PLConLow=ExpERMTooWetG(e)%CondLow(t-1)

                         obj(i)%ecoScorelowExp(e)=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
                         EXIT tloop3
                       END IF  !End if too wet time>m
                     END DO tloop3
                   ELSE !If startScorelowExp(e)<=0.4
             tloop4: DO t=1,ExpERMTooWetP(e)%numlineslow
                       IF (e==4 .and. m<=ExpERMTooWetP(e)%TimeLow(2)) THEN
                         obj(i)%ecoScorelowExp(e)=obj(i-1)%ecoScorelowExp(e)
                       ELSE
                         IF (ExpERMTooWetP(e)%TimeLow(t)>=m) THEN
                           PLInUp=ExpERMTooWetP(e)%TimeLow(t)
                           PLConUp=ExpERMTooWetP(e)%CondLow(t)
                           PLInLow=ExpERMTooWetP(e)%TimeLow(t-1)
                           PLConLow=ExpERMTooWetP(e)%CondLow(t-1)

                           obj(i)%ecoScorelowExp(e)=(REAL(m)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
                           EXIT tloop4
                         END IF  !End if too wet time>m
                       END IF  !End if e==4
                     END DO tloop4
                   END IF !End if prevScore>0.4
                 END IF !End if TooWetCount>1
               END IF ! End if m< too wet start
             END DO  !End from e=1 to numexpertlower

        IF (i==1) THEN
          IF (startScoreup<0.00001) THEN
            obj(i)%ecoScoreup=0.0
          END IF     
          IF (startScorelow<0.00001) THEN
            obj(i)%ecoScorelow=0.0
          END IF  
        END IF 

        IF(i>1) THEN
          IF (obj(i-1)%ecoScoreup<0.00001) THEN
            obj(i)%ecoScoreup=0.0
          END IF
          DO e=1,numexpertupper 
            IF (obj(i-1)%ecoScoreupExp(e)<0.00001) THEN
              obj(i)%ecoScoreupExp(e)=0.0
            END IF
          END DO     
          IF (obj(i-1)%ecoScorelow<0.00001) THEN
            obj(i)%ecoScorelow=0.0
          END IF 
          DO e=1,numexpertlower 
            IF (obj(i-1)%ecoScorelowExp(e)<0.00001) THEN
              obj(i)%ecoScorelowExp(e)=0.0
            END IF
          END DO  
        END IF 
        
     ! WRITE(200,*)obj(i)%ecoScoreup                          !Only use with aggregated upper and lower bounds across experts
     ! WRITE(300,*)obj(i)%ecoScorelow 
      WRITE(210,fmt)(obj(i)%ecoScoreupExp(e),e=1,numexpertupper)
      WRITE(310,fmt)(obj(i)%ecoScorelowExp(e),e=1,numexpertlower)        
        
      i=i+1 
                  
      END DO
     END IF 
   END DO 
! END DO
       
  END SUBROUTINE envFlowProcessor    
  
  
END MODULE processMod

!-------------------------------------------------------------------------------
!PROGRAM
!-------------------------------------------------------------------------------

PROGRAM process
USE processMod
IMPLICIT NONE
  CALL getFlowRainValues
  CALL estGWDepth
  CALL envFlowProcessor
END PROGRAM process

!------------------------------------------------------------------------------------
!FUNCTIONS
! The purpose of these functions are to lookup the specified condition curves using the current 
! condition and duration, given the curves assume a particular starting condition.
!---------------------------------------------------------------------------------- 

FUNCTION DurAdjustDry(LookupTabTime, LookupTabCond,startScore,listnum)
  REAL, INTENT(IN):: LookupTabTime(listnum), LookupTabCond(listnum), startScore
  INTEGER, INTENT(IN)::listnum
  REAL :: PLInUp,PLConUp,PLInLow,PLConLow
  INTEGER :: t

    IF (startScore<1.0) THEN 
      DO t=1,listnum
        IF (LookupTabCond(t)<=startScore) THEN
          PLInUp=LookupTabTime(t)
          PLConUp=LookupTabCond(t)
          PLInLow=LookupTabTime(t-1)
          PLConLow=LookupTabCond(t-1)
            
          DurAdjustDry= PLInLow + ((startScore-PLConLow)*(PLInUp-PLInLow)/(PLConUp-PLConLow))

          EXIT
        END IF
      END DO !end loop through Dry Upper Bound
          
    ELSE 
      DurAdjustDry= 0.0     
    END IF  !End duration adjustment block
  RETURN
END FUNCTION DurAdjustDry

FUNCTION DurAdjustWet(LookupTabTime, LookupTabCond,startScore,listnum)
  REAL, INTENT(IN):: LookupTabTime(listnum), LookupTabCond(listnum), startScore
  INTEGER, INTENT(IN)::listnum
  REAL :: PLInUp,PLConUp,PLInLow,PLConLow
  INTEGER :: t

    IF (startScore>0.25) THEN 
      DO t=1,listnum
        IF (LookupTabCond(t)>=startScore) THEN
          PLInUp=LookupTabTime(t)
          PLConUp=LookupTabCond(t)
          PLInLow=LookupTabTime(t-1)
          PLConLow=LookupTabCond(t-1)
            
          DurAdjustWet= PLInLow + ((startScore-PLConLow)*(PLInUp-PLInLow)/(PLConUp-PLConLow))
 
          EXIT
        END IF
      END DO
          
    ELSE
      PLInUp=LookupTabTime(2)
      PLConUp=LookupTabCond(2)
      PLInLow=LookupTabTime(1)
      PLConLow=LookupTabCond(1)                         
 
      DurAdjustWet= PLInLow - ((PLConLow-startScore)*(PLInUp-PLInLow)/(PLConUp-PLConLow))     
    END IF
  RETURN
END FUNCTION DurAdjustWet

FUNCTION CalcEcoScoreDry(LookupTabTime, LookupTabCond,prevScore,listnum,duration)
  REAL, INTENT(IN):: LookupTabTime(listnum), LookupTabCond(listnum), prevScore
  INTEGER, INTENT(IN)::listnum, duration
  REAL :: PLInUp,PLConUp,PLInLow,PLConLow
  INTEGER :: t
    IF (prevScore<1.0) THEN        
      DO t=1,listnum
        IF (LookupTabTime(t)>=duration) THEN
          PLInUp=LookupTabTime(t)
          PLConUp=LookupTabCond(t)
          PLInLow=LookupTabTime(t-1)
          PLConLow=LookupTabCond(t-1)
            
          CalcEcoScoreDry=(REAL(duration)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
 
          EXIT
        END IF
      END DO  !End Dry Upper Bound block    
    ELSE
      PLInUp=LookupTabTime(2)
      PLConUp=LookupTabCond(2)
      PLInLow=LookupTabTime(1)
      PLConLow=LookupTabCond(1)

      CalcEcoScoreDry=PLConLow-(PLInLow-REAL(duration))*(PLConUp-PLConLow)/(PLInUp-PLInLow)
    END IF  !End if start score<1.0 loop
  RETURN
END FUNCTION CalcEcoScoreDry

FUNCTION CalcEcoScoreWet(LookupTabTime, LookupTabCond,prevScore,listnum,duration)
  REAL, INTENT(IN):: LookupTabTime(listnum), LookupTabCond(listnum), prevScore
  INTEGER, INTENT(IN)::listnum, duration
  REAL :: PLInUp,PLConUp,PLInLow,PLConLow
  INTEGER :: t

    IF (prevScore>0.25) THEN
      DO t=1,listnum
        IF (LookupTabTime(t)>=duration) THEN
          PLInUp=LookupTabTime(t)
          PLConUp=LookupTabCond(t)
          PLInLow=LookupTabTime(t-1)
          PLConLow=LookupTabCond(t-1)
            
          CalcEcoScoreWet=(REAL(duration)-PLInLow)*(PLConUp-PLConLow)/(PLInUp-PLInLow)+PLConLow
 
          EXIT
        END IF
      END DO      
    ELSE
      PLInUp=LookupTabTime(2)
      PLConUp=LookupTabCond(2)
      PLInLow=LookupTabTime(1)
      PLConLow=LookupTabCond(1)

      CalcEcoScoreWet=PLConLow-(PLInLow-REAL(duration))*(PLConUp-PLConLow)/(PLInUp-PLInLow)
    END IF !End if obj(i-1)%ecoScoreup>0.2
  RETURN
END FUNCTION CalcEcoScoreWet
