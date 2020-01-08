## MAGNETOKALMAN ##

## Kinematic model 1 ##

$$ \left \lbrace \array{ u^+ = u \\ s^+ = s + T_su \\ n^+ = n + T_s\xi u \\ \xi^+ = \xi - T_sk(s)u \\ m_s^+ = m_s \\ m_n^+ = m_n \\ m_z^+ = m_z \\ D^+ = D \\ z^+ = z } \right \rbrace $$

## Kinematic model 2 (the one used: yaw rate compensates the curvature of the road) ##

$$ \left \lbrace \array{ u^+ = u \\ s^+ = s + T_su \\ n^+ = n + T_s\xi u \\ \xi^+ = \xi  \\ m_s^+ = m_s \\ m_n^+ = m_n \\ m_z^+ = m_z \\ D^+ = D \\ z^+ = z } \right \rbrace $$

## measure model ##

 $$ \left[ \begin {array}{c} {\it B0x}+1/4\,{\frac {\mu_{0}\, \left(  3.0
\,s \left( t \right) - 3.0\,s_{0}- \left(  \left(  \left| -s \left( t
 \right) +s_{0} \right|  \right) ^{2}+ \left(  \left| -n \left( t
 \right) +n_{0} \right|  \right) ^{2}+ \left(  \left| z \left( t
 \right)  \right|  \right) ^{2} \right) ^{5/2} \left( {\it mx} \left( 
t \right) \cos \left( \theta_{0} \right) -{\it my} \left( t \right) 
\sin \left( \theta_{0} \right) - \left( \sin \left( \theta_{0}
 \right) {\it mx} \left( t \right) +\cos \left( \theta_{0} \right) {
\it my} \left( t \right)  \right) \xi \left( t \right)  \right) 
 \right) }{\pi\, \left(  \left(  \left| -s \left( t \right) +s_{0}
 \right|  \right) ^{2}+ \left(  \left| -n \left( t \right) +n_{0}
 \right|  \right) ^{2}+ \left(  \left| z \left( t \right)  \right| 
 \right) ^{2} \right) ^{5/2}}}+D \left( t \right) {\it B0x}
\\  {\it B0y}+1/4\,{\frac {\mu_{0}\, \left(  3.0\,n
 \left( t \right) - 3.0\,n_{0}- \left(  \left(  \left| -s \left( t
 \right) +s_{0} \right|  \right) ^{2}+ \left(  \left| -n \left( t
 \right) +n_{0} \right|  \right) ^{2}+ \left(  \left| z \left( t
 \right)  \right|  \right) ^{2} \right) ^{5/2} \left( \sin \left( 
\theta_{0} \right) {\it mx} \left( t \right) +\cos \left( \theta_{0}
 \right) {\it my} \left( t \right) + \left( {\it mx} \left( t \right) 
\cos \left( \theta_{0} \right) -{\it my} \left( t \right) \sin \left( 
\theta_{0} \right)  \right) \xi \left( t \right)  \right)  \right) }{
\pi\, \left(  \left(  \left| -s \left( t \right) +s_{0} \right| 
 \right) ^{2}+ \left(  \left| -n \left( t \right) +n_{0} \right| 
 \right) ^{2}+ \left(  \left| z \left( t \right)  \right|  \right) ^{2
} \right) ^{5/2}}}+D \left( t \right) {\it B0y}
\\  {\it B0z}+1/4\,{\frac {\mu_{0}\, \left(  3.0\,z
 \left( t \right) - \left(  \left(  \left| -s \left( t \right) +s_{0}
 \right|  \right) ^{2}+ \left(  \left| -n \left( t \right) +n_{0}
 \right|  \right) ^{2}+ \left(  \left| z \left( t \right)  \right| 
 \right) ^{2} \right) ^{5/2}{\it mz} \left( t \right)  \right) }{\pi\,
 \left(  \left(  \left| -s \left( t \right) +s_{0} \right|  \right) ^{
2}+ \left(  \left| -n \left( t \right) +n_{0} \right|  \right) ^{2}+
 \left(  \left| z \left( t \right)  \right|  \right) ^{2} \right) ^{5/
2}}}+D  {\it B0z}\end {array} \right]  $$

## Considerations ##

* zeros on Q? is that allowed, I expect so since they are affected for other states
* `X(3) = - K*u;    % noise here as yaw rate` Are really using info from K(s) think about this. because the equation is
  $$ \xi  ^+ = \xi + T_s\dot{\psi} - T_sk(s)u $$ and we are using $\dot{\psi}$ as noise.
  I think you are using more $\theta$.
* **Stabilizzazione di corsia: ci vuole qualcosa che contenga l'angolo, tipo barrier function sui lati corsia, o comunque qualcosa che penalizzi pesantemente grosse variazioni di angolo rispetto a quelli "possibili" nella corsia**



## TODO ##

* **controlla sistemi di riferimento dipoli magnetici**
  * ANS: Il modello di misura di ogni magnetometro deve contenere theta (fisso) e xi(t) linearizzato.  **Perciò nel filtro ci devono essere.** Ricordati che però nelle simulazione il theta va messo prima! (o riscrivi funzione)
* **ROTAZIONI SU B_0**
  * Dovrebbero essere corrette!
* controlla funzione simulazione magnetometro (con D) D ~ 1 m^3
  * va

