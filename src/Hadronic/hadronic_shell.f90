module hadronic_shell
    implicit none
    private
    public :: pp_delta
    public :: hic_shell, hic_project
    public :: species_step
    public :: scan_pmax, secondary_rad
    public :: secondary_project
    public :: loss_rates, electron_seq
    public :: photon_loss, effective_time
    public :: pgamma_update
    public :: proton_step, exp_sink
    public :: rate_lum
    public :: project_lum, project_hic
    public :: pair_content
    public :: shell_density
    public :: proc_power, pos_interp, gamma_source
    public :: dist_gev, align_photon, shell_volumes
    public :: shell_geom
contains

! Single-shell pp delta-approximation wrapper.
subroutine pp_delta(Num_p,ep,pden,ntarget, &
                                      Num_gamma,egam,Num_nu,enu,Num_pair,epair, &
                                      kappa,fpion,fneutral,model,qgam, &
                                      qnu,qpair,ploss)
    use hadronic_pp, only: pp_source
    implicit none
    integer, intent(in) :: Num_p,Num_gamma,Num_nu,Num_pair,model
    real(8), intent(in), dimension(Num_p) :: ep,pden
    real(8), intent(in) :: ntarget
    real(8), intent(in), dimension(Num_gamma) :: egam
    real(8), intent(in), dimension(Num_nu) :: enu
    real(8), intent(in), dimension(Num_pair) :: epair
    real(8), intent(in) :: kappa,fpion,fneutral
    real(8), intent(out), dimension(Num_gamma) :: qgam
    real(8), intent(out), dimension(Num_nu) :: qnu
    real(8), intent(out), dimension(Num_pair) :: qpair
    real(8), intent(out), dimension(Num_p) :: ploss

    call pp_source(Num_p,ep,pden,ntarget, &
                   Num_gamma,egam,Num_nu,enu,Num_pair,epair, &
                   model,qgam,qnu,qpair,ploss, &
                   kappa,fpion,fneutral)
end subroutine pp_delta

! Single-shell hadronic inverse-Compton wrapper.
subroutine hic_shell(Num_had,ehad,Num_ph,eph,phden, &
                                         pden,pipden,pimden,mumlden, &
                                         mumrden,muplden,muprden, &
                                         iphmin,epsp,epspi,epsmu, &
                                         cp,cpi,cmu,dln,dep,jmp,depi, &
                                         jmpi,demu,jmmu)
    use hadronic_ic, only: ic_init, &
                                           ic_operator
    implicit none
    integer, intent(in) :: Num_had,Num_ph,iphmin
    real(8), intent(in), dimension(Num_had) :: ehad
    real(8), intent(in), dimension(Num_ph) :: eph,phden
    real(8), intent(in), dimension(Num_had) :: pden,pipden,pimden
    real(8), intent(in), dimension(Num_had) :: mumlden,mumrden
    real(8), intent(in), dimension(Num_had) :: muplden,muprden
    real(8), intent(out), dimension(Num_ph) :: epsp,epspi,epsmu
    real(8), intent(out) :: cp,cpi,cmu,dln
    integer, intent(out), dimension(Num_had) :: dep,jmp,depi,jmpi
    integer, intent(out), dimension(Num_had) :: demu,jmmu

    call ic_init(Num_had,ehad,Num_ph,eph, &
                                                iphmin,dln, &
                                                dep,jmp,depi,jmpi,demu,jmmu)
    call ic_operator(Num_had,Num_ph,phden,pden, &
                                                   pipden,pimden,mumlden, &
                                                   mumrden,muplden, &
                                                   muprden,dln,dep,jmp, &
                                                   depi,jmpi,demu,jmmu,epsp, &
                                                   epspi,epsmu,cp,cpi, &
                                                   cmu)
end subroutine hic_shell

! proton-only hadronic IC：目标光子插值、IC operator 和 lum 投影在 Fortran 内闭合。
subroutine hic_project(num_had,num_ph,num_align,ehad,eph, &
                                     phden,pden,vol,lum)
    implicit none
    integer, intent(in) :: num_had,num_ph,num_align
    real(8), intent(in), dimension(num_had) :: ehad
    real(8), intent(in), dimension(num_ph) :: eph,phden
    real(8), intent(in), dimension(num_had) :: pden
    real(8), intent(in) :: vol
    real(8), intent(out), dimension(num_ph) :: lum
    real(8), dimension(num_align) :: aph,phalign
    real(8), dimension(num_had) :: zero_had
    real(8), dimension(num_align) :: epsilon_p,epsilon_pi,epsilon_mu
    real(8) :: coeff_p,coeff_pi,coeff_mu,dln
    integer, dimension(num_had) :: delta_p,jmp,delta_pi,jmpi
    integer, dimension(num_had) :: delta_mu,jmmu

    zero_had=0d0
    call align_photon(num_had,num_ph,num_align,ehad,eph,aph)
    call pos_interp(num_ph,num_align,eph,phden, &
                                            aph,phalign)
    call hic_shell(num_had,ehad,num_align,aph, &
                                       phalign,pden,zero_had,zero_had,zero_had,zero_had, &
                                       zero_had,zero_had,0,epsilon_p,epsilon_pi,epsilon_mu,coeff_p,coeff_pi, &
                                       coeff_mu,dln,delta_p,jmp,delta_pi,jmpi,delta_mu,jmmu)
    call project_hic(num_align,num_ph,aph,epsilon_p,epsilon_pi,epsilon_mu, &
                                            vol,eph,lum)
end subroutine hic_project

! 单壳层二级强子输运：源项映射、π/μ 推进和 neutron sink 在 Fortran 内闭合。
subroutine species_step(ng,ns,gamma,esrc,qnin, &
                                              qpipin,qpimin, &
                                              qmmlin,qmmrin, &
                                              qmplin,qmprin, &
                                              nlossin,dt,bfield,divrate, &
                                              vol,nprev,pipp,pimp, &
                                              mmlp,mmrp,mplp, &
                                              mprp,nnext,pipn,pimn, &
                                              mmln,mmrn,mpln, &
                                              mprn)
    use constants
    use hadronic_species, only: species_advance
    implicit none
    integer, intent(in) :: ng,ns
    real(8), intent(in), dimension(ng) :: gamma
    real(8), intent(in), dimension(ns) :: esrc
    real(8), intent(in) :: dt,bfield,divrate
    real(8), intent(in), dimension(ns) :: qnin,qpipin
    real(8), intent(in) :: vol
    real(8), intent(in), dimension(ns) :: qpimin,qmmlin
    real(8), intent(in), dimension(ns) :: qmmrin,qmplin
    real(8), intent(in), dimension(ns) :: qmprin,nlossin
    real(8), intent(in), dimension(ng) :: nprev,pipp,pimp
    real(8), intent(in), dimension(ng) :: mmlp,mmrp
    real(8), intent(in), dimension(ng) :: mplp,mprp
    real(8), intent(out), dimension(ng) :: nnext,pipn,pimn
    real(8), intent(out), dimension(ng) :: mmln,mmrn
    real(8), intent(out), dimension(ng) :: mpln,mprn
    integer :: i
    real(8), dimension(ng) :: en,epi,emu
    real(8), dimension(ng) :: qn,qpip,qpim
    real(8), dimension(ng) :: qml,qmr,qpl,qpr
    real(8), dimension(ng) :: nloss,ntran

    do i=1,ng
        en(i)=gamma(i)*Para_m_n_GeV
        epi(i)=gamma(i)*Para_m_pi_charged_GeV
        emu(i)=gamma(i)*Para_m_mu_GeV
    end do
    call gamma_source(ns,ng,esrc,qnin, &
                                      en,Para_m_n_GeV,vol,qn)
    call gamma_source(ns,ng,esrc,qpipin, &
                                      epi,Para_m_pi_charged_GeV,vol,qpip)
    call gamma_source(ns,ng,esrc,qpimin, &
                                      epi,Para_m_pi_charged_GeV,vol,qpim)
    call gamma_source(ns,ng,esrc,qmmlin, &
                                      emu,Para_m_mu_GeV,vol,qml)
    call gamma_source(ns,ng,esrc,qmmrin, &
                                      emu,Para_m_mu_GeV,vol,qmr)
    call gamma_source(ns,ng,esrc,qmplin, &
                                      emu,Para_m_mu_GeV,vol,qpl)
    call gamma_source(ns,ng,esrc,qmprin, &
                                      emu,Para_m_mu_GeV,vol,qpr)
    call species_advance(ng,gamma,dt,bfield,divrate,nprev, &
                                           pipp,pimp,mmlp, &
                                           mmrp,mplp,mprp, &
                                           qn,qpip,qpim,qml, &
                                           qmr,qpl,qpr,ntran, &
                                           pipn,pimn,mmln, &
                                           mmrn,mpln,mprn)
    call pos_interp(ns,ng,esrc,nlossin, &
                                            en,nloss)
    call exp_sink(ng,ntran,nloss,dt,nnext)
end subroutine species_step

! 全局质子最大 Lorentz 因子：逐半径壳层估计并取最大值作为输运网格上界。
subroutine scan_pmax(nr,rad,bulk,bfield,eta_acc,gmaxall)
    use hadronic_accel, only: gamma_limit
    implicit none
    integer, intent(in) :: nr
    real(8), intent(in), dimension(nr) :: rad,bulk,bfield
    real(8), intent(in) :: eta_acc
    real(8), intent(out) :: gmaxall
    integer :: i
    real(8), dimension(2) :: gscan,extrate
    real(8) :: gmax,gdyn,gsyn,gext
    logical :: has_extlim

    if (eta_acc <= 0d0) error stop "hadronic eta_acc must be positive."
    gscan=[1d0,2d0]
    extrate=[1d0,2d0]
    gmaxall=0d0
    do i=1,nr
        call gamma_limit("proton",bfield(i),rad(i),bulk(i),eta_acc, &
                                         2,gscan,extrate,.false.,gmax,gdyn, &
                                         gsyn,gext,has_extlim)
        gmaxall=max(gmaxall,gmax)
    end do
    if (gmaxall <= 1d0) error stop "hadronic maximum proton Lorentz factor must exceed unity."
end subroutine scan_pmax

! Single-shell secondary-radiation wrapper.
subroutine secondary_rad(Num_had,ehad,Num_ph,eph,pipden, &
                                                 pimden,mumlden,mumrden, &
                                                 muplden,muprden,phden, &
                                                 iphmin,bmag,pisynr, &
                                                 musynr,piicr,muicr)
    use hadronic_secondary, only: secondary_calc
    implicit none
    integer, intent(in) :: Num_had,Num_ph,iphmin
    real(8), intent(in), dimension(Num_had) :: ehad
    real(8), intent(in), dimension(Num_ph) :: eph
    real(8), intent(in), dimension(Num_had) :: pipden,pimden
    real(8), intent(in), dimension(Num_had) :: mumlden,mumrden
    real(8), intent(in), dimension(Num_had) :: muplden,muprden
    real(8), intent(in), dimension(Num_ph) :: phden
    real(8), intent(in) :: bmag
    real(8), intent(out), dimension(Num_ph) :: pisynr,musynr
    real(8), intent(out), dimension(Num_ph) :: piicr,muicr
    real(8) :: dsyn,dic
    real(8), dimension(Num_ph,Num_had) :: kpion,kmuon
    integer, dimension(Num_had) :: depi,jmpi,demu,jmmu

    call secondary_calc(Num_had,ehad,Num_ph,eph,pipden,pimden,mumlden,mumrden, &
                        muplden,muprden,phden,iphmin,bmag, &
                        pisynr,musynr, &
                        piicr,muicr,dsyn, &
                        kpion,kmuon,dic,depi,jmpi,demu,jmmu)
end subroutine secondary_rad

! 投影后的二级 π/μ 辐射：分布映射、目标光子插值、辐射和 lum 投影在 Fortran 内闭合。
subroutine secondary_project(num_had,num_ph,num_align,ehad, &
                                                     eph,phden, &
                                                     pippg,pimpg, &
                                                     mmlpg, &
                                                     mmrpg, &
                                                     mplpg, &
                                                     mprpg, &
                                                     vol,bmag,pisynlum, &
                                                     musynlum,piiclum,muiclum)
    use constants
    implicit none
    integer, intent(in) :: num_had,num_ph,num_align
    real(8), intent(in), dimension(num_had) :: ehad
    real(8), intent(in), dimension(num_ph) :: eph,phden
    real(8), intent(in), dimension(num_had) :: pippg,pimpg
    real(8), intent(in), dimension(num_had) :: mmlpg,mmrpg
    real(8), intent(in), dimension(num_had) :: mplpg,mprpg
    real(8), intent(in) :: vol,bmag
    real(8), intent(out), dimension(num_ph) :: pisynlum,musynlum
    real(8), intent(out), dimension(num_ph) :: piiclum,muiclum
    integer :: i
    real(8), dimension(num_had) :: gsp,epi,emu
    real(8), dimension(num_align) :: aph
    real(8), dimension(num_align) :: phalign
    real(8), dimension(num_had) :: pipden,pimden
    real(8), dimension(num_had) :: mmlgev,mmrgev,mplgev,mprgev
    real(8), dimension(num_align) :: pisyn,musyn,piic,muic

    do i=1,num_had
        gsp(i)=ehad(i)/Para_m_p_GeV
        epi(i)=gsp(i)*Para_m_pi_charged_GeV
        emu(i)=gsp(i)*Para_m_mu_GeV
    end do
    call dist_gev(num_had,num_had,epi,pippg, &
                                          ehad,Para_m_pi_charged_GeV,vol,pipden)
    call dist_gev(num_had,num_had,epi,pimpg, &
                                          ehad,Para_m_pi_charged_GeV,vol,pimden)
    call dist_gev(num_had,num_had,emu,mmlpg, &
                                          ehad,Para_m_mu_GeV,vol,mmlgev)
    call dist_gev(num_had,num_had,emu,mmrpg, &
                                          ehad,Para_m_mu_GeV,vol,mmrgev)
    call dist_gev(num_had,num_had,emu,mplpg, &
                                          ehad,Para_m_mu_GeV,vol,mplgev)
    call dist_gev(num_had,num_had,emu,mprpg, &
                                          ehad,Para_m_mu_GeV,vol,mprgev)
    call align_photon(num_had,num_ph,num_align,ehad,eph,aph)
    call pos_interp(num_ph,num_align,eph,phden, &
                                            aph,phalign)
    call secondary_rad(num_had,ehad,num_align,aph,pipden, &
                                               pimden,mmlgev,mmrgev,mplgev, &
                                               mprgev,phalign,0,bmag, &
                                               pisyn,musyn,piic,muic)
    call project_lum(num_align,num_ph,aph,pisyn, &
                                                  vol,eph,pisynlum)
    call project_lum(num_align,num_ph,aph,musyn, &
                                                  vol,eph,musynlum)
    call project_lum(num_align,num_ph,aph,piic, &
                                                  vol,eph,piiclum)
    call project_lum(num_align,num_ph,aph,muic, &
                                                  vol,eph,muiclum)
end subroutine secondary_project

! 连续冷却率 wrapper：绝热项加同步项，可选 quantum synch 修正。
subroutine loss_rates(ng,gamma,bfield,tdyn,mass,quantum_syn,ltot)
    use constants
    use hadronic_base, only: quant_factor
    implicit none
    integer, intent(in) :: ng,quantum_syn
    real(8), intent(in), dimension(ng) :: gamma
    real(8), intent(in) :: bfield,tdyn,mass
    real(8), intent(out), dimension(ng) :: ltot
    integer :: i
    real(8) :: csyn,mratio,syn_loss

    if (bfield < 0d0) error stop "hadronic continuous loss rates require bfield >= 0."
    if (tdyn <= 0d0) error stop "hadronic continuous loss rates require tdyn > 0."
    mratio=mass/Para_m_e_GeV
    csyn=Para_sigmaT*bfield*bfield/(6d0*pi*Para_m_e*Para_c)/(mratio**3)
    do i=1,ng
        syn_loss=csyn*gamma(i)*gamma(i)
        if (quantum_syn /= 0) syn_loss=syn_loss*quant_factor(gamma(i),bfield,mass)
        ltot(i)=gamma(i)/tdyn+syn_loss
    end do
end subroutine loss_rates

! BH/pp 二级电子序列：按半径壳层推进冷却谱，并输出同步辐射源项。
subroutine electron_seq(num_e,num_nu,nr,gamma_e,rad,bulk,bfield, &
                                                   freq,src,synidx,nthr,quantum_syn, &
                                                   eden,lsyn,seed,srcrad)
    use constants
    use electron_radiation_kernel, only: syn_state
    use hadronic_remap, only: remap_loggamma
    implicit none
    integer, intent(in) :: num_e,num_nu,nr,synidx,nthr,quantum_syn
    real(8), intent(in), dimension(num_e) :: gamma_e
    real(8), intent(in), dimension(nr) :: rad,bulk,bfield
    real(8), intent(in), dimension(num_nu) :: freq
    real(8), intent(in), dimension(num_e,nr) :: src
    real(8), intent(out), dimension(num_e,nr) :: eden
    real(8), intent(out), dimension(num_nu,nr) :: lsyn,seed
    real(8), intent(out), dimension(num_e,nr) :: srcrad
    integer :: ir
    real(8), dimension(num_e) :: ltot,prev,next
    real(8) :: dr,dt,tdyn
    real(8), dimension(num_nu) :: ptmp,tautmp

    eden=0d0; lsyn=0d0; seed=0d0; srcrad=0d0
    prev=0d0
    do ir=1,nr
        call shell_geom(nr,rad,bulk,ir,dr,dt)
        tdyn=rad(ir)/(bulk(ir)*Para_c)
        call loss_rates(num_e,gamma_e,bfield(ir),tdyn, &
                                               Para_m_e_GeV,quantum_syn,ltot)
        if (ir == 1) then
            next=prev
            srcrad(:,ir)=0d0
        else
            call remap_loggamma(num_e,gamma_e,prev, &
                                                        src(:,ir)*dt,ltot,dt,next)
            srcrad(:,ir)=src(:,ir)*dt/dr*gamma_e(:)
        end if
        call syn_state(synidx,rad(ir),bfield(ir),num_e,num_nu,nthr, &
                                    gamma_e,next,freq,ptmp,lsyn(:,ir), &
                                    seed(:,ir),tautmp)
        eden(:,ir)=next
        prev=next
    end do
end subroutine electron_seq

! 光子损失闭合：由壳层 comoving path time 得到 tau，并返回局部 survival 平均因子。
subroutine photon_loss(num_ph,nr,rad,bulk,ish,rate,tau, &
                                           survival)
    use constants
    implicit none
    integer, intent(in) :: num_ph,nr,ish
    real(8), intent(in), dimension(nr) :: rad,bulk
    real(8), intent(in), dimension(num_ph) :: rate
    real(8), intent(out), dimension(num_ph) :: tau,survival
    integer :: iph
    real(8) :: dr,dt

    if (any(rate < 0d0)) error stop "hadronic photon loss closure requires non-negative loss rate."
    call shell_geom(nr,rad,bulk,ish,dr,dt)
    do iph=1,num_ph
        tau(iph)=rate(iph)*dt
        if (tau(iph) > 1d-6) then
            survival(iph)=(1d0-dexp(-tau(iph)))/tau(iph)
        else if (tau(iph) > 0d0) then
            survival(iph)=1d0-tau(iph)/2d0+tau(iph)*tau(iph)/6d0
        else
            survival(iph)=1d0
        end if
    end do
end subroutine photon_loss

! 相互作用有效时间：对指数 sink 的同一步 reinjection 积分。
subroutine effective_time(nrate,rate,dt,teff)
    use constants
    implicit none
    integer, intent(in) :: nrate
    real(8), intent(in), dimension(nrate) :: rate
    real(8), intent(in) :: dt
    real(8), intent(out), dimension(nrate) :: teff
    integer :: i
    real(8) :: tau

    if (dt < 0d0) error stop "hadronic interaction effective time requires dt >= 0."
    if (any(rate < 0d0)) error stop "hadronic interaction effective time requires non-negative rates."
    do i=1,nrate
        if (rate(i) > 0d0) then
            tau=rate(i)*dt
            if (tau < 1d-4) then
                teff(i)=dt*(1d0-tau/2d0+tau*tau/6d0)
            else
                teff(i)=(1d0-dexp(-tau))/rate(i)
            end if
        else
            teff(i)=dt
        end if
    end do
end subroutine effective_time

! pγ 后质子更新：指数 sink 与同一步 reinjection 在同一壳层内闭合。
subroutine pgamma_update(ng,dtran,loss,reinj, &
                                            vol,dt,next)
    use constants
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: dtran,loss,reinj
    real(8), intent(in) :: vol,dt
    real(8), intent(out), dimension(ng) :: next
    integer :: i
    real(8), dimension(ng) :: efftime

    if (vol <= 0d0) error stop "hadronic p-gamma proton update requires positive shell volume."
    call effective_time(ng,loss,dt,efftime)
    do i=1,ng
        next(i)=dtran(i)*dexp(-loss(i)*dt)+ &
                   efftime(i)*vol*reinj(i)*Para_m_p_GeV
    end do
end subroutine pgamma_update

! 单壳层质子输运：连续冷却 + BH/pp 损失 + pγ sink/reinjection 在 Fortran 内闭合。
subroutine proton_step(ng,gamma,prev,qinj,bfield,tdyn,mass, &
                                             quantum_syn,bhloss,pploss,pgloss,pgreinj, &
                                             vol,dt,next)
    use hadronic_remap, only: remap_loggamma
    implicit none
    integer, intent(in) :: ng,quantum_syn
    real(8), intent(in), dimension(ng) :: gamma,prev,qinj
    real(8), intent(in) :: bfield,tdyn,mass,vol,dt
    real(8), intent(in), dimension(ng) :: bhloss,pploss,pgloss
    real(8), intent(in), dimension(ng) :: pgreinj
    real(8), intent(out), dimension(ng) :: next
    real(8), dimension(ng) :: ltot,dtran

    call loss_rates(ng,gamma,bfield,tdyn,mass,quantum_syn,ltot)
    ltot(1:ng)=ltot(1:ng)+bhloss(1:ng)+pploss(1:ng)
    call remap_loggamma(ng,gamma,prev,qinj,ltot,dt,dtran)
    call pgamma_update(ng,dtran,pgloss,pgreinj, &
                                          vol,dt,next)
end subroutine proton_step

! 指数 sink：用于单步粒子损失项 N -> N exp(-rate dt)。
subroutine exp_sink(nval,values,loss,dt,vnext)
    use constants
    implicit none
    integer, intent(in) :: nval
    real(8), intent(in), dimension(nval) :: values,loss
    real(8), intent(in) :: dt
    real(8), intent(out), dimension(nval) :: vnext
    integer :: i

    if (dt < 0d0) error stop "hadronic exponential sink requires dt >= 0."
    if (any(loss < 0d0)) error stop "hadronic exponential sink requires non-negative rates."
    do i=1,nval
        vnext(i)=values(i)*dexp(-loss(i)*dt)
    end do
end subroutine exp_sink

! 将每能量反应率谱转换为 shell lum 谱。
subroutine rate_lum(ne,en,rate,vol, &
                                                   lum)
    use constants
    implicit none
    integer, intent(in) :: ne
    real(8), intent(in), dimension(ne) :: en,rate
    real(8), intent(in) :: vol
    real(8), intent(out), dimension(ne) :: lum

    if (vol <= 0d0) error stop "hadronic lum conversion requires positive shell volume."
    lum(1:ne)=vol*rate(1:ne)* &
                             en(1:ne)*Para_h_GeV*Para_GeV2erg
end subroutine rate_lum

! rate 谱投影：先转壳层 lum，再映射到目标 photon energy grid。
subroutine project_lum(ns,nd,esrc,rsrc, &
                                                    vol,edst,ldst)
    use constants
    implicit none
    integer, intent(in) :: ns,nd
    real(8), intent(in), dimension(ns) :: esrc,rsrc
    real(8), intent(in), dimension(nd) :: edst
    real(8), intent(in) :: vol
    real(8), intent(out), dimension(nd) :: ldst
    real(8), dimension(ns) :: lsrc

    if (vol <= 0d0) error stop "hadronic lum projection requires positive shell volume."
    lsrc(1:ns)=vol*rsrc(1:ns)*esrc(1:ns)* &
                              Para_h_GeV*Para_GeV2erg
    call pos_interp(ns,nd,esrc,lsrc, &
                                            edst,ldst)
end subroutine project_lum

! hadronic IC lum 投影：合并 p/pi/mu IC 能量源项并映射到目标 photon grid。
subroutine project_hic(ns,nd,esrc,epsp,epspi, &
                                              epsmu,vol,edst,ldst)
    use constants
    implicit none
    integer, intent(in) :: ns,nd
    real(8), intent(in), dimension(ns) :: esrc,epsp,epspi,epsmu
    real(8), intent(in), dimension(nd) :: edst
    real(8), intent(in) :: vol
    real(8), intent(out), dimension(nd) :: ldst
    real(8), dimension(ns) :: lsrc

    if (vol <= 0d0) error stop "hadronic IC lum projection requires positive shell volume."
    lsrc(1:ns)=vol*(epsp(1:ns)+epspi(1:ns)+ &
                              epsmu(1:ns))*Para_h_GeV*Para_GeV2erg
    call pos_interp(ns,nd,esrc,lsrc, &
                                            edst,ldst)
end subroutine project_hic

! BH/pp 二级电子源项：把每 GeV 产生率合并为壳层内每 gamma 注入率。
! Convert pp rates plus both BH lepton charges per GeV into shell-integrated rates per gamma.
subroutine pair_content(num_e,ppair,bhpair,include_pp, &
                                           include_bh,vol,src)
    use constants
    implicit none
    integer, intent(in) :: num_e,include_pp,include_bh
    real(8), intent(in), dimension(num_e) :: ppair,bhpair
    real(8), intent(in) :: vol
    real(8), intent(out), dimension(num_e) :: src

    if (vol <= 0d0) error stop "hadronic pair source content requires positive shell volume."
    src=0d0
    if (include_pp /= 0) src(1:num_e)=src(1:num_e)+ppair(1:num_e)
    if (include_bh /= 0) src(1:num_e)=src(1:num_e)+2d0*bhpair(1:num_e)
    src(1:num_e)=vol*Para_m_e_GeV*src(1:num_e)
end subroutine pair_content

! 壳层粒子密度归一化：从每 gamma 数量转换为每 GeV 体密度。
subroutine shell_density(ng,dgamma,mass,vol, &
                                             dgev)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: dgamma
    real(8), intent(in) :: mass,vol
    real(8), intent(out), dimension(ng) :: dgev

    if (mass <= 0d0) error stop "hadronic shell density per GeV requires positive particle mass."
    if (vol <= 0d0) error stop "hadronic shell density per GeV requires positive shell volume."
    dgev(1:ng)=dgamma(1:ng)/(vol*mass)
end subroutine shell_density

! AM3 分过程功率归并：积分每个过程 lum，并按质子能量分布投到 hadron grid。
subroutine proc_power(num_had,nproce,num_process,ehad,dn_had, &
                                     eproc,rproc,vol,ppow)
    use constants
    implicit none
    integer, intent(in) :: num_had,nproce,num_process
    real(8), intent(in), dimension(num_had) :: ehad,dn_had
    real(8), intent(in), dimension(nproce) :: eproc
    real(8), intent(in), dimension(num_process,nproce) :: rproc
    real(8), intent(in) :: vol
    real(8), intent(out), dimension(num_process,num_had) :: ppow
    integer :: i,j
    real(8), dimension(num_had) :: pweight
    real(8), dimension(nproce) :: lum
    real(8) :: wtot,ptot

    if (vol <= 0d0) error stop "hadronic process power requires positive shell volume."
    do i=1,num_had
        pweight(i)=dn_had(i)*ehad(i)
    end do
    wtot=trapz(num_had,ehad,pweight)
    ppow=0d0
    if (wtot <= 0d0) return
    do j=1,num_process
        lum(1:nproce)=vol*rproc(j,1:nproce)* &
                                      eproc(1:nproce)*Para_h_GeV*Para_GeV2erg
        ptot=trapz(nproce,eproc,lum)
        ppow(j,1:num_had)=pweight(1:num_had)/wtot*ptot
    end do
contains
    real(8) function trapz(num_x,x,y)
        implicit none
        integer, intent(in) :: num_x
        real(8), intent(in), dimension(num_x) :: x,y
        integer :: k

        trapz=0d0
        do k=1,num_x-1
            trapz=trapz+0.5d0*(y(k)+y(k+1))*(x(k+1)-x(k))
        end do
    end function trapz
end subroutine proc_power

! 正值 log-log 插值：只使用 finite 且正的源点，范围外输出零。
subroutine pos_interp(ns,nd,xs,ys,xd,yd)
    use ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer, intent(in) :: ns,nd
    real(8), intent(in), dimension(ns) :: xs,ys
    real(8), intent(in), dimension(nd) :: xd
    real(8), intent(out), dimension(nd) :: yd
    integer :: i,j,nvalid
    real(8), dimension(ns) :: xv,yv
    real(8) :: frac

    yd=0d0; nvalid=0
    do i=1,ns
        if (ieee_is_finite(xs(i)) .and. ieee_is_finite(ys(i)) .and. xs(i) > 0d0 .and. ys(i) > 0d0) then
            nvalid=nvalid+1
            xv(nvalid)=xs(i)
            yv(nvalid)=ys(i)
        end if
    end do
    if (nvalid < 2) return
    do i=1,nd
        if (xd(i) < xv(1) .or. xd(i) > xv(nvalid) .or. .not. ieee_is_finite(xd(i))) cycle
        do j=1,nvalid-1
            if (xd(i) >= xv(j) .and. xd(i) <= xv(j+1)) then
                frac=dlog(xd(i)/xv(j))/dlog(xv(j+1)/xv(j))
                yd(i)=dexp(dlog(yv(j))+frac*dlog(yv(j+1)/yv(j)))
                exit
            end if
        end do
    end do
end subroutine pos_interp

! pγ 源项映射：从每 GeV 反应率谱变换到二级粒子的每 gamma 注入率。
subroutine gamma_source(ns,nd,esrc,qsrc,edst, &
                                        mass,vol,qdst)
    implicit none
    integer, intent(in) :: ns,nd
    real(8), intent(in), dimension(ns) :: esrc,qsrc
    real(8), intent(in), dimension(nd) :: edst
    real(8), intent(in) :: mass,vol
    real(8), intent(out), dimension(nd) :: qdst

    if (mass <= 0d0) error stop "hadronic source per gamma requires positive particle mass."
    if (vol <= 0d0) error stop "hadronic source per gamma requires positive shell volume."
    call pos_interp(ns,nd,esrc,qsrc, &
                                            edst,qdst)
    qdst(1:nd)=vol*mass*qdst(1:nd)
end subroutine gamma_source

! 二级粒子分布映射：从壳层内每 gamma 数量变换为局部每 GeV 数密度谱。
subroutine dist_gev(ns,nd,esrc,dsrc,edst, &
                                            mass,vol,ddst)
    implicit none
    integer, intent(in) :: ns,nd
    real(8), intent(in), dimension(ns) :: esrc,dsrc
    real(8), intent(in), dimension(nd) :: edst
    real(8), intent(in) :: mass,vol
    real(8), intent(out), dimension(nd) :: ddst

    if (mass <= 0d0) error stop "hadronic distribution per GeV requires positive particle mass."
    if (vol <= 0d0) error stop "hadronic distribution per GeV requires positive shell volume."
    call pos_interp(ns,nd,esrc,dsrc, &
                                            edst,ddst)
    ddst(1:nd)=ddst(1:nd)/(vol*mass)
end subroutine dist_gev

! 按 hadron log spacing 构造对齐 photon energy grid。
subroutine align_photon(num_had,num_ph,num_out,ehad,eph, &
                                           aphgev)
    implicit none
    integer, intent(in) :: num_had,num_ph,num_out
    real(8), intent(in), dimension(num_had) :: ehad
    real(8), intent(in), dimension(num_ph) :: eph
    real(8), intent(out), dimension(num_out) :: aphgev
    integer :: i
    real(8) :: dlnhad,lmin

    if (num_had < 2 .or. num_ph < 2 .or. num_out < 1) error stop "hadronic aligned photon grid needs valid grids."
    dlnhad=dlog(ehad(2)/ehad(1))
    lmin=dlog(eph(1))
    do i=1,num_out
        aphgev(i)=dexp(lmin+dlnhad*dble(i-1))
    end do
end subroutine align_photon

! 相邻半径间对 dt'/dR=1/(beta Gamma c) 作梯形积分；首点是零时长初态。
! Trapezoid-integrate proper time between radii; the first point is the zero-duration initial state.
subroutine shell_geom(nr,rad,bulk,ir,dr,dt)
    use constants
    implicit none
    integer, intent(in) :: nr,ir
    real(8), intent(in), dimension(nr) :: rad,bulk
    real(8), intent(out) :: dr,dt
    real(8) :: u_prev,u_now

    if (ir == 1) then
        dr=0d0
        dt=0d0
        return
    end if
    if (bulk(ir-1) <= 1d0 .or. bulk(ir) <= 1d0) &
        error stop "hadronic sequence shell dt requires bulk > 1."
    dr=rad(ir)-rad(ir-1)
    if (dr <= 0d0) error stop "hadronic sequence shell radii must be strictly increasing."
    u_prev=dsqrt((bulk(ir-1)-1d0)*(bulk(ir-1)+1d0))
    u_now=dsqrt((bulk(ir)-1d0)*(bulk(ir)+1d0))
    dt=0.5d0*dr/Para_c*(1d0/u_prev+1d0/u_now)
end subroutine shell_geom

! 半径壳层体积：第一个壳层以内边界 R=0，之后使用相邻半径。
subroutine shell_volumes(nr,rad,vol)
    use constants
    implicit none
    integer, intent(in) :: nr
    real(8), intent(in), dimension(nr) :: rad
    real(8), intent(out), dimension(nr) :: vol
    integer :: i
    real(8) :: rprev

    if (nr < 1) error stop "hadronic shell volumes require at least 1d0 radius."
    if (rad(1) <= 0d0) error stop "hadronic shell radii must be positive."
    vol(1)=(4d0/3d0)*pi*rad(1)**3
    rprev=rad(1)
    do i=2,nr
        if (rad(i) <= 0d0) error stop "hadronic shell radii must be positive."
        if (rad(i) <= rad(i-1)) error stop "hadronic shell radii must be strictly increasing."
        vol(i)=(4d0/3d0)*pi*(rad(i)**3-rprev**3)
        rprev=rad(i)
    end do
end subroutine shell_volumes

end module hadronic_shell
