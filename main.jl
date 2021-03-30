import Pkg;
Pkg.activate(@__DIR__);
using LinearAlgebra, StaticArrays
using RigidBodyDynamics, RigidBodySim
using MeshCat, MeshCatMechanisms
vis = Visualizer();open(vis)
using Gadfly, Cairo, Fontconfig

q0 = [0.01;-0.5;-0.0;-2.0;-0.3;1.5;-0.7;0.1;0.1]
qd = [0.0;0.0;0.0;0.0;0.0;pi;0.01;0.01;0.01]

function display_urdf(urdfPath,vis)
    # Displays mechanism at config all zeros
    # urdfPath must be a string
    mechanism = parse_urdf(Float64,urdfPath)
    state = MechanismState(mechanism)
    zero_velocity!(state)
    set_configuration!(state,q0)
    mvis = MechanismVisualizer(mechanism, URDFVisuals(urdfPath),vis)
#    manipulate!(state) do x
#        set_configuration!(mvis, configuration(x))
#    end
#    for bd in bodies(mechanism)
#        setelement!(mvis,default_frame(bd),0.5,"$bd")
#    end
    return mvis, mechanism
end

display_urdf("panda.urdf",vis)

# Assignment 3 ##############################################


function traj(t::Float64)
    # compute the desired joint angle at time t

    # Setup
    A = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 2 0 0 0;
         1 t t^2 t^3 t^4 t^5;
         0 1 2*t 3*t^2 4*t^3 5*t^4;
         0 0 2 6*t 12*t^2 20*t^3]

    b1 = [q0[1];0;0;qd[1];0;0]
    b2 = [q0[2];0;0;qd[2];0;0]
    b3 = [q0[3];0;0;qd[3];0;0]
    b4 = [q0[4];0;0;qd[4];0;0]
    b5 = [q0[5];0;0;qd[5];0;0]
    b6 = [q0[6];0;0;qd[6];0;0]
    b7 = [q0[7];0;0;qd[7];0;0]
    b8 = [q0[8];0;0;qd[8];0;0]
    b9 = [q0[9];0;0;qd[9];0;0]

    # solving for constants of polynomials

    x1 = A\b1
    x2 = A\b2
    x3 = A\b3
    x4 = A\b4
    x5 = A\b5
    x6 = A\b6
    x7 = A\b7
    x8 = A\b8
    x9 = A\b9

    # constructing qf, vf, and af at time t

    qf1 = [1 t t^2 t^3 t^4 t^5]*x1
    qf2 = [1 t t^2 t^3 t^4 t^5]*x2
    qf3 = [1 t t^2 t^3 t^4 t^5]*x3
    qf4 = [1 t t^2 t^3 t^4 t^5]*x4
    qf5 = [1 t t^2 t^3 t^4 t^5]*x5
    qf6 = [1 t t^2 t^3 t^4 t^5]*x6
    qf7 = [1 t t^2 t^3 t^4 t^5]*x7
    qf8 = [1 t t^2 t^3 t^4 t^5]*x8
    qf9 = [1 t t^2 t^3 t^4 t^5]*x9

    qdotf1 = [0 1 2*t 3*t^2 4*t^3 5*t^4]*x1
    qdotf2 = [0 1 2*t 3*t^2 4*t^3 5*t^4]*x2
    qdotf3 = [0 1 2*t 3*t^2 4*t^3 5*t^4]*x3
    qdotf4 = [0 1 2*t 3*t^2 4*t^3 5*t^4]*x4
    qdotf5 = [0 1 2*t 3*t^2 4*t^3 5*t^4]*x5
    qdotf6 = [0 1 2*t 3*t^2 4*t^3 5*t^4]*x6
    qdotf7 = [0 1 2*t 3*t^2 4*t^3 5*t^4]*x7
    qdotf8 = [0 1 2*t 3*t^2 4*t^3 5*t^4]*x8
    qdotf9 = [0 1 2*t 3*t^2 4*t^3 5*t^4]*x9

    qddotf1 = [0 0 2 6*t 12*t^2 20*t^3]*x1
    qddotf2 = [0 0 2 6*t 12*t^2 20*t^3]*x2
    qddotf3 = [0 0 2 6*t 12*t^2 20*t^3]*x3
    qddotf4 = [0 0 2 6*t 12*t^2 20*t^3]*x4
    qddotf5 = [0 0 2 6*t 12*t^2 20*t^3]*x5
    qddotf6 = [0 0 2 6*t 12*t^2 20*t^3]*x6
    qddotf7 = [0 0 2 6*t 12*t^2 20*t^3]*x7
    qddotf8 = [0 0 2 6*t 12*t^2 20*t^3]*x8
    qddotf9 = [0 0 2 6*t 12*t^2 20*t^3]*x9

    qf = [qf1;qf2;qf3;qf4;qf5;qf6;qf7;qf8;qf9]

    qdotf = [qdotf1;qdotf2;qdotf3;qdotf4;qdotf5;qdotf6;qdotf7;qdotf8;qdotf9]

    qddotf = [qddotf1;qddotf2;qddotf3;qddotf4;qddotf5;qddotf6;qddotf7;qddotf8;qddotf9]

    return qf, qdotf, qddotf
end


function control_PD!(τ, t, state)
    # Compute a value for τ
    τ .= -diagm([100,80, 80, 50, 20, 20, 20,10,30]) * velocity(state) - diagm([150,150, 75, 75, 50, 50, 30,30,25])*(configuration(state) - [0.0;0.0;0.0;0.0;0.0;pi;0.01;0.01;0.01])
    # Saturate
    act_sat = 50; # Actuator limits
    τ .= map( x -> x > act_sat ? act_sat : x,τ)
    τ .= map( x -> x < -act_sat ? -act_sat : x,τ)
end

function control_CTC!(τ, t, state)
    # Compute a value for τ

    τ .= mass_matrix*similar(velocity(state)) + dynamics_bias

    # Saturate
    act_sat = 50; # Actuator limits
    τ .= map( x -> x > act_sat ? act_sat : x,τ)
    τ .= map( x -> x < -act_sat ? -act_sat : x,τ)
end


qf, qdotf, qddotf = traj(10.)

# set panda to initial configuration

display_urdf("panda.urdf",vis)

sol1 = simulate(traj(10.),10.,control_PD!)

# return panda to initial configuration

display_urdf("panda.urdf",vis)

sol2 = simulate(traj(10.),10.,control_CTC!)
