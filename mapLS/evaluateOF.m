function fitnessvalue=evaluateOF(connectivity)

    smallN = 1e-12;

    MaxDispRef = 1;
    lambdaMass = 10;

    omega=100;
    plotvar=0;
    MAXmass=2e-6;
    Dmin=50e-6;
    % E,nu,rho
    beamMaterialProperties=[1e12,0.49,1000];
    % E,nu,rho,thickness
    membProperties=[10e10,.49,1000,2e-6];

    try
        flclear fem

        % COMSOL version
        clear vrsn
        vrsn.name = 'COMSOL 3.2';
        vrsn.ext = 'a';
        vrsn.major = 0;
        vrsn.build = 300;
        vrsn.rcs = '$Name:  $';
        vrsn.date = '$Date: 2005/12/20 19:02:30 $';
        fem.version = vrsn;

        % create Geometry
        %[g3D lengthB]=creategeometry(nTrees,delta0,rF0,dF0,Y0,strTree0);
        [geom1,g3D]=createGeometry(connectivity);

        % Geometry 1
        f.objs={g3D};
        f.name={'menbrane'};
        f.tags={'mb'};


        fem.draw=struct('f',f);
        fem.geom=geomcsg(fem);

        % Initialize mesh for geometry 1
        fem.mesh=meshinit(fem, ...
                          'hmaxfact',0.35, ...
                          'hcurve',0.3, ...
                          'hgrad',1.35, ...
                          'hcutoff',0.005, ...
                          'hnarrow',0.85);

        % (Default values are not included)

        % Application mode 1

        % appl
        clear appl
        appl.mode.class = 'Sme3DEulerBeam';
        appl.module = 'SME';
        appl.gporder = 6;
        appl.cporder = 1;
        appl.assignsuffix = '_smeul3d';

        % prop
        clear prop
        prop.analysis='freq';
        appl.prop = prop;

        % edg
        clear edg
        edg.heightz = 'diam';
        edg.heighty = 'diam';

        edg.A = 'pi*diam^2/4';
        edg.Iyy = 'pi*(diam/2)^4/4';
        edg.Izz = 'pi*(diam/2)^4/4';
        edg.J = 'pi*(diam/2)^4/2';

        edg.Fz = 'Pressure*diam';

        edg.localxp = -.001;
        edg.localyp = .001;

        edg.rho = 'mat2_rho';
        edg.E = 'mat2_E';
        edg.nu = 'mat2_nu';

        NoEdges = flgeomnes(fem.geom);
        edg.ind = ones(1,NoEdges);
        appl.edg = edg;

            clear pnt
        pnt.Hthy = {0,1};
        pnt.Hz = {0,1};
        pnt.Hthz = {0,1};
        pnt.Hy = {0,1};
        pnt.Hx = {0,1};
        pnt.Hthx = {0,1};
        pnt.ind = ones(1,flgeomnv(fem.geom));
        vert=flgeomvtx(fem.geom);
        aux1=find(vert(1,:)==0);
        for i=1:length(aux1)
            pnt.ind(aux1(i))=2;
        end
        appl.pnt = pnt;

        % update fem with Euler beam settings
        fem.appl{1} = appl;

        % Application mode 2

        % appl
        clear appl
        appl.mode.class = 'SmeShell';
        appl.module = 'SME';
        appl.gporder = 2;
        appl.cporder = 1;
        appl.assignsuffix = '_smsh';

        % prop
        clear prop
        prop.analysis='freq';
        appl.prop = prop;

        % edg
        clear edg
        edg.Hthy = {0,1};
        edg.Hz = {0,1};
        edg.Hthz = {0,1};
        edg.Hy = {0,1};
        edg.Hx = {0,1};
        edg.Hthx = {0,1};
        edg.ind = ones(1,NoEdges);
        edgeG = flgeomse(fem.geom);
        vtxG = flgeomvtx(fem.geom);

        % Setting Constraints
        for iEd=1:NoEdges
             vtxE = edgeG(:,iEd);
             xEd0 = vtxG(1,vtxE(1));
             xEd1 = vtxG(1,vtxE(2));
             if(xEd0<smallN && xEd1<smallN)
                 iEdType = 2;
             else
                 iEdType = 1;
             end;
             edg.ind(iEd) = iEdType;
        end;
        appl.edg = edg;

        % bnd
        clear bnd
        bnd.rho = 'mat1_rho';
        bnd.Fz = 'Pressure';
        bnd.E = 'mat1_E';
        bnd.nu = 'mat1_nu';
        bnd.thickness = 'mat1_thickness';
        bnd.ind = ones(1,flgeomnbs(fem.geom));
        appl.bnd = bnd;

        % update fem with sheel settings
        fem.appl{2} = appl;
        fem.frame = {'ref'};
        fem.border = 1;
        fem.units = 'SI';

        % Global expressions
        fem.globalexpr = {'Pressure','-100000',...
            'diam','0.005'};

        % Library materials

        clear lib

        % membrane
        lib.mat{1}.name='Membrane';
        lib.mat{1}.varname='mat1';
        lib.mat{1}.variables.rho='1000';
        lib.mat{1}.variables.nu='0.49';
        lib.mat{1}.variables.E='10e10';
        lib.mat{1}.variables.thickness='2e-6';

        % veins
        lib.mat{2}.name='Veins';
        lib.mat{2}.varname='mat2';
        lib.mat{2}.variables.rho='1000';
        lib.mat{2}.variables.nu='0.49';
        lib.mat{2}.variables.E='1e12';


        fem.lib = lib;

        % Multiphysics
        fem=multiphysics(fem);

        % Extend mesh
        fem.xmesh=meshextend(fem);

        % Solve problem
        fem.sol=femlin(fem, ...
                    'symmetric','off', ...
                    'solcomp',{'thx','thz','w','u','thy','v'}, ...
                    'outcomp',{'thx','thz','w','u','thy','v'}, ...
                    'thresh',0.1,...
                    'uscale','init',...
                    'linsolver','umfpack');

        save femstr fem;

        if (plotvar==1)
            postplot(fem, ...
                 'tridata',{'w','cont','internal'}, ...
                 'trimap','jet', ...
                 'maxminbnd','w', ...
                 'deformscale',1, ...
                 'deformsub',{'','',''}, ...
                 'grid','on', ...
                 'campos',[-0.070913581287553,-0.0828832318453904,0.0668530716820182], ...
                 'camtarget',[6.32746417950751E-5,0.00961560768623077,-4.61486902565408E-4], ...
                 'camup',[0,0,1], ...
                 'camva',11.437140256824206);
        end

        maxdisp=max([abs(postmax(fem,'w','edim',1)),abs(postmax(fem,'w','edim',2)),...
                     abs(postmin(fem,'w','edim',1)),abs(postmin(fem,'w','edim',2))]);

        mass=postint(fem,'rho_smeul3d*A_smeul3d','edim',1);

        fitnessvalue=maxdisp/MaxDispRef+0*lambdaMass*max([(mass-MAXmass)/MAXmass,0]);
    catch
        fitnessvalue=1e6;
    end
    
end


function [geom1,geom2]=createGeometry(C)
    % Geometry 2
    k=1;
    g{k}=rect2(0.01,0.01,'base','corner','pos',[0,0]);
    matrix = C.matrix;
    coord = C.coord;
    Nvert = C.last;
    for i=1:Nvert
        aux = [zeros(1,i),matrix(i,i+1:end)];
        col = find(aux==1);
        for j=1:length(col);
            k=k+1;
            carr={curve2([coord(i,1),coord(col(j),1)],[coord(i,2),coord(col(j),2)],[1,1])};
            g{k}=geomcoerce('curve',carr);
        end
    end        
    
    geom1=geomcoerce('solid',g);
    geom2=embed(geom1,'Wrkpln',[0 1 0;0 0 1;0 0 0]);  

end