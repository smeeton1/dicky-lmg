
def IntPhi(b,a,omega,omega0,eta,gamma,i):
  phi=(b-a)*(numpy.arccos(math.sqrt((omega*omega0)*((eta/(2*omega0))*pow(b,2)+b-i)/(2*pow(gamma,2)*(1-pow(b,2)))))+numpy.arccos(math.sqrt((omega*omega0)*((eta/(2*omega0))*pow(a,2)+a-i)/(2*pow(gamma,2)*(1-pow(a,2))))))/2
  return phi


def DoSR1(xmin,xmax,e1,e2,omega,omega0,eta,gamma):
  CDoS = []
  XDoS= []
  i=e1+0.01
  XDoS.append(i)
  while i<e2:
    z1=-(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))+omega*omega0)/(4*pow(gamma,2)+eta*omega)
    z2=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    i+=0.1
    XDoS.append(i)
    phi=IntPhi(z2,z1,omega,omega0,eta,gamma,i)
    CDoS.append((z2+1)/2+phi/math.pi)
  while e2<=i<xmax
    CDoS.append(1)
    i+=0.1
    XDoS.append(i)
  return CDoS,XDoS

def DoSR2(xmin,xmax,e1,e2,emin,omega,omega0,eta,gamma):
  CDoS = []
  XDoS= []
  i=e1+0.01
  while i<emin:
    zn=-(math.sqrt(omega0*(omega-2*i*eta))+omega)/eta
    zp=(math.sqrt(omega0*(omega-2*i*eta))-omega)/eta
    i+=0.1
    XDoS.append(i)
    phi=IntPhi(zn,zp,omega,omega0,eta,gamma,i)
    CDoS.append(phi/math.pi)
  while emin<=i<e2
    z2=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    zp=(math.sqrt(omega0*(omega-2*i*eta))-omega)/eta
    phi=IntPhi(z2,zp,omega,omega0,eta,gamma,i)
    CDoS.append((z2+1)/2+phi/math.pi)
    i+=0.1
    XDoS.append(i)
  while e2<=i<xmax
    CDoS.append(1)
    i+=0.1
    XDoS.append(i)
  return CDoS,XDoS


def DoSR3(xmin,xmax,e1,e2,emin,eNe,omega,omega0,eta,gamma):
  CDoS = []
  XDoS= []
  i=e1+0.01
  XDoS.append(i)
  while i<eNe:
    zn=-(math.sqrt(omega0*(omega-2*i*eta))+omega)/eta
    zp=(math.sqrt(omega0*(omega-2*i*eta))-omega)/eta
    i+=0.1
    XDoS.append(i)
    phi=IntPhi(zn,zp,omega,omega0,eta,gamma,i)
    CDoS.append(phi/math.pi)
  while emNe<=i<emin
    z1=-(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))+omega*omega0)/(4*pow(gamma,2)+eta*omega)
    z2=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    zp=(math.sqrt(omega0*(omega-2*i*eta))-omega)/eta
    zn=-(math.sqrt(omega0*(omega-2*i*eta))+omega)/eta
    phi=IntPhi(zn,z1,omega,omega0,eta,gamma,i)
    phi2=IntPhi(z2,zp,omega,omega0,eta,gamma,i)
    CDoS.append((z2-z1)/2+(phi+phi2)/math.pi)
    i+=0.1
    XDoS.append(i)
  while emin<=i<e2
    z2=(math.sqrt(16*pow(gamma,4)+4*pow(gamma,2)*omega*(eta+2*i*omega0)+omega0*pow(omega,2)*(2*i*eta+omega0))-omega*omega0)/(4*pow(gamma,2)+eta*omega)
    zp=(math.sqrt(omega0*(omega-2*i*eta))-omega)/eta
    phi=IntPhi(z2,zp,omega,omega0,eta,gamma,i)
    CDoS.append((z2+1)/2+phi/math.pi)
    i+=0.1
    XDoS.append(i)
  while e2<=i<xmax
    CDoS.append(1)
    i+=0.1
    XDoS.append(i)
  return CDoS,XDoS



