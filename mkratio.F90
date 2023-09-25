! #undef OPENMP
#define OPENMP
PROGRAM mkratio

!=======================================================
! PFT LAI ratio of tree/grass
!=======================================================

   USE netcdf
   USE omp_lib

   IMPLICIT NONE

   INTEGER, parameter :: r8 = selected_real_kind(12)
   INTEGER, parameter :: xydim = 1200
   INTEGER, parameter :: nday = 46
   INTEGER, parameter :: lon_points=86400
   INTEGER, parameter :: lat_points=43200
   INTEGER, parameter :: iday = 24 ! LAI数据, 第iday个8天
   
   CHARACTER (len=255) :: LAI_DIR = "/tera06/yuanhua/modis/global_lai_15s/" ! V6.1
   CHARACTER (len=255) :: VCF_DIR = "/tera06/yuanhua/mksrf/vcf_15s/" ! V6.1
   CHARACTER (len=255) :: REGFILE = "reg_5x5test"
   CHARACTER (len=4)   :: year    = "2001"
   CHARACTER (len=255) :: filename

   ! short data
   INTEGER          , dimension(4)  :: reg
   CHARACTER (len=4), dimension(4)  :: creg
   INTEGER :: siglb,eiglb,sjglb,ejglb, i_omp,num_threads
   INTEGER :: iglb,jglb, j,i, j1,i1, win, cnt
   
   REAL :: lai,grasum,tresum,passpct
   REAL :: fillvalue, gracnt,trecnt, pctt, pcth, pctb !,summ
   
   INTEGER  :: XY2D(2)
   INTEGER  :: argn
   INTEGER  :: t1,t2,t3,t4,t5
   LOGICAL  :: flag

   ! input vars
   REAL, dimension(86400, 43200) :: pcttdata, pcthdata, summ !, pctbdata
   REAL, dimension(86400, 43200) :: laidata
  
   REAL, dimension(43200) :: latdata
   REAL, dimension(86400) :: londata

   ! input vars id
   INTEGER :: ncid, pcttid, pcthid, laiid !, pctbid

   ! ***阈值设定***
   ! 网格内树、草覆盖比例的阈值: 树或草的比例大于90%
   ! 从当前像素开始搜索, 满足条件的格子数超过10个,则停止向外搜索;窗口最大不超过maxwin
!    INTEGER , parameter    :: maxwin   = 241  ! 1200/240 = 5°/1°
   INTEGER , parameter    :: maxwin   = 721
   INTEGER , parameter    :: fracthre = 90   ! the threshold of tree or grass fraction
   REAL    , parameter    :: ratothre = 10   ! the threshold of the number of grids that satisfy fracthre


   ! output data  
   REAL, dimension(1200)      :: lats, lons
   REAL, dimension(1200,1200) :: laio!, pctto, pctho
   REAL, allocatable, dimension(:,:) :: tregrarat

   ! output vars id
   INTEGER :: lat_dimid, lat_vid, lon_dimid, lon_vid, tregrarat_id, laio_id!, pctto_id, pctho_id


   allocate( tregrarat(1200,1200))
   tregrarat(:,:) = 0.
   call system_clock(t1)
   write(*,*) 't1: ', t1, 'ms'

   ! get args from command line
   argn = IARGC()
   IF (argn > 0) THEN
    CALL getarg(1, REGFILE)
    CALL getarg(2, year)
   ENDIF


   PRINT*, '>>> Read raw data...'
   ! 读 VCF=====================================================
   filename = TRIM(VCF_DIR)//'vcf-modis-'//TRIM(year)//'.nc'
   PRINT*, filename

   CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )

   CALL check( nf90_inq_varid(ncid, 'PCTT', pcttid  ) )
   CALL check( nf90_get_var  (ncid, pcttid, pcttdata) )

   CALL check( nf90_inq_varid(ncid, 'PCTH', pcthid  ) )
   CALL check( nf90_get_var  (ncid, pcthid, pcthdata) )

   ! CALL check( nf90_inq_varid(ncid, 'PCTB', pctbid  ) )
   ! CALL check( nf90_get_var  (ncid, pctbid, pctbdata) )

   CALL check( nf90_inq_varid(ncid, 'lat', lat_vid  ) )
   CALL check( nf90_get_var  (ncid, lat_vid, latdata) )

   CALL check( nf90_inq_varid(ncid, 'lon', lon_vid  ) )
   CALL check( nf90_get_var  (ncid, lon_vid, londata) )

   CALL check( nf90_close(ncid) )


! !    TEST------------------------------------------------
!    jglb = 67746
!    iglb = 14480
!    DO j=jglb-2,jglb+2
!       DO i=iglb-2,iglb+2
!          print*, 'j=',j,'i=',i, 'pctt=',pcttdata(j,i), 'pcth=',pcthdata(j,i)
!       ENDDO
!    ENDDO
! !    TEST------------------------------------------------

   ! VCF前处理
   pcttdata = pcttdata/0.8_r8
   pcthdata = pcthdata/0.8_r8
   summ = pcttdata + pcthdata
   where ( summ > 100 )
      pcttdata = pcttdata * 100./summ
      pcthdata = pcthdata * 100./summ
   end where
   PRINT*, 'global VCF preprocessed completed'   
   call system_clock(t2)
   write(*,*) 't2: ', t2, 'ms'
   write(*,*) 'time usage in preprocessing vcf:', (t2-t1)/(1000.*60.), 'min'

   ! 读 LAI=====================================================
   filename = TRIM(LAI_DIR)//'global_lai_15s_'//TRIM(year)//'.nc'
   PRINT*, filename
   CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )
   CALL check( nf90_inq_varid(ncid, 'lai' , laiid       ) )
   CALL check( nf90_get_var  (ncid, laiid , laidata(:,:), &
                    start=(/1,1,iday/)                  , &
                    count=(/lon_points,lat_points,1/)   ) )
   CALL check( nf90_close(ncid) )
!    print*, 'test lai',laidata(86400,43200)

   ! assign
   laidata (:,:) = laidata (:,:)*0.1_r8

   PRINT*, 'raw data read completed'   
   call system_clock(t3)
   write(*,*) 't3: ', t3, 'ms'
   write(*,*) 'time usage in reading data:', (t3-t2)/(1000.*60.), 'min'

   ! Reg5x5 loop=================================================
   open(13,FILE=trim(REGFILE))
   open(14,FILE=trim(REGFILE))

   DO WHILE(.TRUE.)
      
      READ(13,*,END=200) reg
      READ(14,*,END=201) creg
      PRINT*, 'region:', creg

      For/ output
      lats(:) = latdata( (90 - reg(1))  / 5 * 1200 + 1 : (90 - reg(1))  / 5 * 1200 + 1200 )
      lons(:) = londata( (reg(2) + 180) / 5 * 1200 + 1 : (reg(2) + 180) / 5 * 1200 + 1200 )
      laio(:,:) = laidata( (reg(2) + 180) / 5 * 1200 + 1 : (reg(2) + 180) / 5 * 1200 + 1200, &
         (90 - reg(1))  / 5 * 1200 + 1 : (90 - reg(1))  / 5 * 1200 + 1200)
!       pctto(:,:) = pcttdata( (reg(2) + 180) / 5 * 1200 + 1 : (reg(2) + 180) / 5 * 1200 + 1200, &
!          (90 - reg(1))  / 5 * 1200 + 1 : (90 - reg(1))  / 5 * 1200 + 1200)
!       pctho(:,:) = pcthdata( (reg(2) + 180) / 5 * 1200 + 1 : (reg(2) + 180) / 5 * 1200 + 1200, &
!          (90 - reg(1))  / 5 * 1200 + 1 : (90 - reg(1))  / 5 * 1200 + 1200)
   
      tregrarat(:,:) = -1.      
      passpct = 0

#ifdef OPENMP 
      ! loop for each small 500m grid
      call omp_set_num_threads(92)
      num_threads = omp_get_max_threads()
      print *, "线程数：", num_threads

!$OMP PARALLEL DO &! NUM_THREADS(92) SCHEDULE(STATIC, 1) &
!$OMP PRIVATE(i,j,iglb,jglb) &
!$OMP PRIVATE(i1,j1,win) &
!$OMP PRIVATE(trecnt,gracnt,tresum,grasum,sjglb,ejglb,siglb,eiglb) &
!$OMP PRIVATE(pctt,pcth,pctb,lai,summ,flag,cnt)
#endif

      DO j=1,xydim,1
         DO i=1,xydim,1
!       DO j=546,546,1
!          DO i=80,80,1
            flag = .false.

            ! 把i,j映射到全球            
            iglb = (90 - reg(1))  / 5 * 1200 + i
            jglb = (reg(2) + 180) / 5 * 1200 + j
            
            IF ( laidata(jglb,iglb) > 1e-6 ) THEN  ! 该网格的LAI不是缺省值

               !先在最大窗口中判断一下
               IF ( count(pcttdata(jglb-(maxwin-1)/2 : jglb+(maxwin-1)/2, &
                                   iglb-(maxwin-1)/2 : iglb+(maxwin-1)/2) .GT. fracthre) &
                         .LT. ratothre ) CYCLE

               IF ( count(pcthdata(jglb-(maxwin-1)/2 : jglb+(maxwin-1)/2, &
                                   iglb-(maxwin-1)/2 : iglb+(maxwin-1)/2) .GT. fracthre) &
                         .LT. ratothre ) CYCLE
               
! #ifdef OPENMP 
!                i_omp = omp_get_thread_num()
!                IF (i_omp == 48) THEN
!                   print*, 'Thread', i_omp, 'is working'
!                   print*, i, j
!                ENDIF
! #endif              
               win = 1    ! 如果自己满足阈值要求的话, 也算在内
               trecnt = 0
               gracnt = 0  
               tresum = 0
               grasum = 0  
               DO WHILE(.not. flag)  ! 确定窗口的边缘
                  ! 把窗口内的数据遍历             
                  sjglb = jglb-(win-1)/2
                  ejglb = jglb+(win-1)/2
                  siglb = iglb-(win-1)/2
                  eiglb = iglb+(win-1)/2   
                  
                  cnt = 0 
                  j1 = sjglb
                  i1 = siglb
                  DO WHILE( cnt .lt. (win**2 - (win-2)**2) )  ! 窗口边边的循环
                     pctt = pcttdata(j1,i1)
                     pcth = pcthdata(j1,i1)
                     ! pctb = 100. - summ
                     
                     lai  = laidata (j1,i1) 

                     ! 判断是否满足条件，对满足条件的格子进行计数
                     IF (pctt .GT. fracthre) THEN
                        trecnt = trecnt + 1
                        tresum = tresum + lai
                     ENDIF

                     IF (pcth .GT. fracthre) THEN
                        gracnt = gracnt + 1
                        grasum = grasum + lai
                     ENDIF

                     cnt = cnt + 1
                     IF (j1 .lt. ejglb .and. i1 .eq. siglb) THEN
                        j1 = j1 + 1
                     ELSE IF (i1 .lt. eiglb .and. j1 .eq. ejglb) THEN
                        i1 = i1 + 1
                     ELSE IF (j1 .gt. sjglb .and. i1 .eq. eiglb) THEN
                        j1 = j1 - 1
                     ELSE IF (i1 .gt. siglb .and. j1 .eq. sjglb) THEN
                        i1 = i1 - 1                     
                     ENDIF

                  ENDDO

                  ! 判断满足条件的格子数是否超过预定义的阈值
                  IF ( (trecnt .GT. ratothre) .AND. (gracnt .GT. ratothre) ) THEN
                     ! 计算树草比例并退出
                     tregrarat(j,i) = (tresum/trecnt) / (grasum/gracnt)
!                      print*, 'pctt=',pcttdata(jglb-(win-1)/2:jglb+(win-1)/2,iglb-(win-1)/2:iglb+(win-1)/2)
!                      print*, 'pcth=',pcthdata(jglb-(win-1)/2:jglb+(win-1)/2,iglb-(win-1)/2:iglb+(win-1)/2)
!                      print*, 'lai=',laidata(jglb-(win-1)/2:jglb+(win-1)/2,iglb-(win-1)/2:iglb+(win-1)/2)
!                      print*, 'j=',j, 'i=',i, 'pctt=',pctt, 'pcth=',pcth
!                      print*, 'tresum=',tresum,'trecnt=',trecnt,'grasum=',grasum,'gracnt=',gracnt
!                      print*, 'j=',j, 'i=',i, 'lairatio=',tregrarat(j,i),'win=',win,'frac=',fracthre
                     passpct = passpct + 1
                     flag = .true.
                  ELSE
                     ! 扩大窗口
                     win = win + 2
                  ENDIF
                  
                  IF (win > maxwin) THEN
                     flag = .true.
                  ENDIF

               ENDDO

            ENDIF
         ENDDO
      ENDDO

#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      call system_clock(t4)
      write(*,*) 't4: ', t4, 'ms'
      write(*,*) 'time usage for calculation:', (t4-t3)/(1000.*60.), 'min'
      
      print*, 'Percent of available ratio = ', passpct/1200./1200.*100

#ifdef OPENMP
      filename = './ratio/ompRG_'//TRIM(creg(1))//'_'//&
                  TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.MOD'//TRIM(year)//'.nc'
#else
      filename = './ratio/sinRG_'//TRIM(creg(1))//'_'//&
                  TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.MOD'//TRIM(year)//'.nc'
#endif
      
      PRINT*, filename
      CALL check( nf90_create(filename, NF90_NETCDF4, ncid) )

      ! define dimensions
      CALL check( nf90_def_dim(ncid, "lat", xydim, lat_dimid) )
      CALL check( nf90_def_dim(ncid, "lon", xydim, lon_dimid) )

      ! define&put variables
      ! ---------------------------------------
      XY2D = (/lon_dimid, lat_dimid/)
      fillvalue = -1.
      CALL check( nf90_def_var(ncid, "lat"          , NF90_FLOAT, lat_dimid, lat_vid, deflate_level=6) )
      CALL check( nf90_def_var(ncid, "lon"          , NF90_FLOAT, lon_dimid, lon_vid, deflate_level=6) )
      CALL check( nf90_def_var(ncid, "LAI"          , NF90_FLOAT,  XY2D, laio_id, deflate_level=6) )
      CALL check( nf90_def_var(ncid, "TRE_GRA_RATIO", NF90_DOUBLE, XY2D, tregrarat_id, deflate_level=6))

      CALL check( nf90_put_att(ncid, lat_vid, "long_name", "Latitude"     ) )
      CALL check( nf90_put_att(ncid, lat_vid, "units"    , "degrees_north") )
      CALL check( nf90_put_att(ncid, lon_vid, "long_name", "Longitude"    ) )
      CALL check( nf90_put_att(ncid, lon_vid, "units"    , "degrees_east" ) )
      CALL check( nf90_put_att(ncid, laio_id             , "units"     , "-"    ) )
      CALL check( nf90_put_att(ncid, laio_id             , "long_name" , "Leaf area index") )
      CALL check( nf90_put_att(ncid, tregrarat_id        , "units"     , "-"    ) )
      CALL check( nf90_put_att(ncid, tregrarat_id        , "long_name" , "LAI ratio of tree/grass") )
      ! NetCDF: Not a valid data type or _FillValue type mismatch
      ! CALL check( nf90_put_att(ncid, tregrarat_id          , "_FillValue", fillvalue    ) )

      CALL check( nf90_enddef(ncid) )

      CALL check( nf90_inq_varid(ncid, "lat", lat_vid ) )
      CALL check( nf90_put_var  (ncid, lat_vid, lats  ) )
      CALL check( nf90_inq_varid(ncid, "lon", lon_vid ) )
      CALL check( nf90_put_var  (ncid, lon_vid, lons  ) )
      CALL check( nf90_inq_varid(ncid, "LAI", laio_id ) )
      CALL check( nf90_put_var  (ncid, laio_id   , laio  ) )
      CALL check( nf90_inq_varid(ncid, "TRE_GRA_RATIO", tregrarat_id ) )
      CALL check( nf90_put_var  (ncid, tregrarat_id   , tregrarat  ) )
      CALL check( nf90_close(ncid) )

   ENDDO

   200 close(13)
   201 close(14)

   call system_clock(t5)
   write(*,*) 't5: ', t5, 'ms'
   write(*,*) 'time usage for create nc:', (t5-t4)/(1000.*60.), 'min'
   write(*,*) 'time usage for all:', (t5-t1)/(1000.*60.), 'min'


CONTAINS

   SUBROUTINE check(status)
      INTEGER, intent(in) :: status

      IF (status /= nf90_noerr) THEN
         print *, trim( nf90_strerror(status))
         stop 2
      ENDIF
   END SUBROUTINE check

END PROGRAM mkratio
