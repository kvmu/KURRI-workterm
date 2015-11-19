![alt tag](https://www.rri.kyoto-u.ac.jp/wp-content/themes/rri/image/etoptitle.jpg)

# Kyoto University Research Reactor Institute (KURRI)
Kyoto University Research Reactor Institute (KURRI) was established in 1963 for the joint use program among Japanese universities to promote the research and education in the fields of nuclear energy and radiation application. Two nuclear reactors, the Kyoto University research Reactor (KUR) and the Kyoto University Critical Assembly (KUCA), and related research facilities have been being used since then, and nowadays greater expectations are being put on the research and education activities at our institute for the issues of energy and environment and for the innovative applications of radiation [[homepage]].

# Overview
This repository is a collection of some of the work I did during my time at KURRI as a Research Assistant. My supervisors were Dr. Yoshiharu Mori (Mori-sensei) [[Google Scholar]], and Dr. Yoshihiro Ishi (Ishi-sensei). I worked in the Innovation Laboratory and worked with the 150-MeV FFAG (Fixed Field Alternating Gradient) Accelerator, a part of the KUCA.  

Below, I provide a small outline of each of the folders in this repository. 

#####Any headers marked with an astersisk (*) are important.

_Special thanks to all the students: Kiwamu Ishii, Yuya Horita, Ryo Masaki, and Abe Yuma for making my stay enjoyable and fun._
_Special thanks to my supervisors and Uesegi-san, for answering all my questions._
_Special thanks to Ms. Yuka Ono for being very organized and always reminding me to pay my rent._

### full_simulation *
A work in progress that is meant to be a transverse-tracking code including the effects of beam-matter interactions. The idea is to be able to simulate using _real_ magnetic fields generated from TOSCA or Pre-Processor, the effect of the carbon-stripping foil to the transverse emittances. 

Currently, it's still in the translation from C++ (Author is Yuya Horita) to Python -- then I plan to add the framework so that it takes into account the effect of the carbon foil.

### long_dynamics_simulation
A simple longitudinal dynamics simulaton and
plots:

1. the phi-W phasespace according to a Hamiltonian (can see on plot)
2. and the longitudinal tracking of a few particles in the phi-W phasespace -- visualizing synchrotron oscillation.

The purpose of this code was for me to understand the background and build some intuition of longitudinal dynamics in accelerators.

### resources
This folder contains some of the pdfs I read/skimmed in order to learn more about beam physics and accelerator physics in preparation for working at KURRI.

### side_project_kiwamu
A fun little comparison of FORTRAN Euler method
code vs Python RK4 code on a system of ODEs. The situation is a parent
isotope is being bombarded by negative muons, which then are captured. The parent splits into daughter
isotopes + neutrons. The application is to use this process to take radioactive
isotopes with a very long half - life to daughter isotopes that are more
stable / less radioactive; safe nuclear waste disposal.

### width_delay_beam_experiment *
This folder contains a full analysis code for a set of beam monitor date (voltage-time series). A report is also included in this folder to provide more details of the calculation.

License
----

MIT


[homepage]: <http://www.rri.kyoto-u.ac.jp/en/outline/preface>
[Google Scholar]: <https://scholar.google.com/citations?user=emRPCooAAAAJ&hl=en>
