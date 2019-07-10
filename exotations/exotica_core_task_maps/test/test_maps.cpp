//
// Copyright (c) 2018, University of Edinburgh
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of  nor the names of its contributors may be used to
//    endorse or promote products derived from this software without specific
//    prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#include <gtest/gtest.h>

#include <exotica_core/exotica_core.h>

// TODO(#437): Activate once solution for pointer casting/dynamic loading is found.
// #include <exotica_core_task_maps/JointAccelerationBackwardDifference.h>
// #include <exotica_core_task_maps/JointJerkBackwardDifference.h>
// #include <exotica_core_task_maps/JointVelocityBackwardDifference.h>

// Extend testing printout //////////////////////

namespace testing
{
namespace internal
{
enum GTestColor
{
    COLOR_DEFAULT,
    COLOR_RED,
    COLOR_GREEN,
    COLOR_YELLOW
};

extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
}
}
#define PRINTF(...)                                                                        \
    do                                                                                     \
    {                                                                                      \
        testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); \
        testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__);    \
    } while (0)

// C++ stream interface
class TestCout : public std::stringstream
{
public:
    ~TestCout()
    {
        PRINTF("%s\n", str().c_str());
    }
};

#define TEST_COUT TestCout()

//////////////////////////////////////////////

using namespace exotica;

std::string urdf_string_ = "<robot name=\"test_robot\"><link name=\"base\"><visual><geometry><cylinder length=\"0.3\" radius=\"0.2\"/></geometry><origin xyz=\"0 0 0.15\"/></visual></link><link name=\"link1\"> <inertial><mass value=\"0.2\"/><origin xyz=\"0 0 0.1\"/><inertia ixx=\"0.00381666666667\" ixy=\"0\" ixz=\"0\" iyy=\"0.0036\" iyz=\"0\" izz=\"0.00381666666667\"/></inertial><visual><geometry><cylinder length=\"0.15\" radius=\"0.05\"/></geometry><origin xyz=\"0 0 0.075\"/></visual> </link><link name=\"link2\"><inertial><mass value=\"0.2\"/><origin xyz=\"0 0 0.1\"/><inertia ixx=\"0.00381666666667\" ixy=\"0\" ixz=\"0\" iyy=\"0.0036\" iyz=\"0\" izz=\"0.00381666666667\"/></inertial><visual><geometry><cylinder length=\"0.35\" radius=\"0.05\"/></geometry><origin xyz=\"0 0 0.175\"/></visual> </link><link name=\"link3\"><inertial><mass value=\"0.2\"/><origin xyz=\"0 0 0.1\"/><inertia ixx=\"0.00381666666667\" ixy=\"0\" ixz=\"0\" iyy=\"0.0036\" iyz=\"0\" izz=\"0.00381666666667\"/></inertial><visual><geometry><cylinder length=\"0.45\" radius=\"0.05\"/></geometry><origin xyz=\"0 0 0.225\"/></visual></link><link name=\"endeff\"><inertial><mass value=\"0.2\"/><origin xyz=\"0 0 0.1\"/><inertia ixx=\"0.00381666666667\" ixy=\"0\" ixz=\"0\" iyy=\"0.0036\" iyz=\"0\" izz=\"0.00381666666667\"/></inertial><visual><geometry><cylinder length=\"0.05\" radius=\"0.1\"/></geometry><origin xyz=\"0 0 -0.025\"/></visual></link><joint name=\"joint1\" type=\"revolute\"><parent link=\"base\"/><child link=\"link1\"/><origin xyz=\"0 0 0.3\" rpy=\"0 0 0\" /><axis xyz=\"0 0 1\" /><limit effort=\"200\"  velocity=\"1.0\" lower=\"-0.4\" upper=\"0.4\"/><safety_controller k_position=\"30\" k_velocity=\"30\" soft_lower_limit=\"-0.4\" soft_upper_limit=\"0.4\"/></joint><joint name=\"joint2\" type=\"revolute\"><parent link=\"link1\"/><child link=\"link2\"/><origin xyz=\"0 0 0.15\" rpy=\"0 0 0\" /><axis xyz=\"0 1 0\" /><limit effort=\"200\"  velocity=\"1.0\" lower=\"-0.4\" upper=\"0.4\"/><safety_controller k_position=\"30\" k_velocity=\"30\" soft_lower_limit=\"-0.4\" soft_upper_limit=\"0.4\"/></joint><joint name=\"joint3\" type=\"revolute\"><parent link=\"link2\"/><child link=\"link3\"/><origin xyz=\"0 0 0.35\" rpy=\"0 0 0\" /><axis xyz=\"0 1 0\" /><limit effort=\"200\"  velocity=\"1.0\" lower=\"-0.4\" upper=\"0.4\"/><safety_controller k_position=\"30\" k_velocity=\"30\" soft_lower_limit=\"-0.4\" soft_upper_limit=\"0.4\"/></joint><joint name=\"joint4\" type=\"fixed\"><parent link=\"link3\"/><child link=\"endeff\"/><origin xyz=\"0 0 0.45\" rpy=\"0 0 0\" /></joint></robot>";
std::string srdf_string_ = "<robot name=\"test_robot\"><group name=\"arm\"><chain base_link=\"base\" tip_link=\"endeff\" /></group><virtual_joint name=\"world_joint\" type=\"fixed\" parent_frame=\"world_frame\" child_link=\"base\" /><group_state name=\"zero\" group=\"arm\"><joint name=\"joint1\" value=\"0\" /><joint name=\"joint2\" value=\"0.3\" /><joint name=\"joint3\" value=\"0.55\" /></group_state><disable_collisions link1=\"base\" link2=\"link1\" reason=\"Adjacent\" /><disable_collisions link1=\"endeff\" link2=\"link3\" reason=\"Adjacent\" /><disable_collisions link1=\"link1\" link2=\"link2\" reason=\"Adjacent\" /><disable_collisions link1=\"link2\" link2=\"link3\" reason=\"Adjacent\" /></robot>";

constexpr bool print_debug_information_ = false;
constexpr int num_trials_ = 100;

bool test_random(UnconstrainedEndPoseProblemPtr problem)
{
    Eigen::VectorXd x(3);
    TEST_COUT << "Testing random configurations:";
    for (int i = 0; i < num_trials_; ++i)
    {
        x = problem->GetScene()->GetKinematicTree().GetRandomControlledState();
        problem->Update(x);
        if (print_debug_information_)
        {
            TEST_COUT << "x = " << x.transpose();
            TEST_COUT << "y = " << problem->Phi.data.transpose();
            TEST_COUT << "jacobian = \n"
                      << problem->jacobian;
        }
    }
    return true;
}

bool test_random(UnconstrainedTimeIndexedProblemPtr problem)
{
    TEST_COUT << "Testing random configurations:";
    for (int i = 0; i < num_trials_; ++i)
    {
        for (int t = 0; t < problem->GetT(); ++t)
        {
            problem->Update(problem->GetScene()->GetKinematicTree().GetRandomControlledState(), t);
        }
    }
    return true;
}

bool test_values(Eigen::MatrixXdRefConst Xref, Eigen::MatrixXdRefConst Yref, Eigen::MatrixXdRefConst Jref, UnconstrainedEndPoseProblemPtr problem, const double eps = 1e-5)
{
    TEST_COUT << "Testing set points:";
    int N = Xref.cols();
    int M = Yref.cols();
    int L = Xref.rows();
    for (int i = 0; i < L; ++i)
    {
        Eigen::VectorXd x = Xref.row(i);
        TaskSpaceVector y = problem->cost.y;
        y.data = Yref.row(i);
        Eigen::MatrixXd jacobian = Jref.middleRows(i * M, M);
        problem->Update(x);
        double errY = (y - problem->Phi).norm();
        double errJ = (jacobian - problem->jacobian).norm();
        if (errY > eps)
        {
            TEST_COUT << "y:  " << problem->Phi.data.transpose();
            TEST_COUT << "y*: " << y.data.transpose();
            ADD_FAILURE() << "Task space error out of bounds: " << errY;
        }
        if (errJ > eps)
        {
            TEST_COUT << "x: " << x.transpose();
            TEST_COUT << "J*:\n"
                      << jacobian;
            TEST_COUT << "J:\n"
                      << problem->jacobian;
            ADD_FAILURE() << "Jacobian error out of bounds: " << errJ;
        }
    }

    return true;
}

bool test_jacobian(UnconstrainedEndPoseProblemPtr problem, const double eps = 1e-5)
{
    constexpr double h = 1e-5;  // NB: Not, this differs from the h for the time-indexed Jacobian

    TEST_COUT << "Testing Jacobian with h=" << h << ", eps=" << eps;
    for (int j = 0; j < num_trials_; ++j)
    {
        Eigen::VectorXd x0(problem->N);
        x0 = problem->GetScene()->GetKinematicTree().GetRandomControlledState();
        problem->Update(x0);
        const TaskSpaceVector y0(problem->Phi);
        const Eigen::MatrixXd J0(problem->jacobian);
        Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(J0.rows(), J0.cols());
        for (int i = 0; i < problem->N; ++i)
        {
            Eigen::VectorXd x(x0);
            x(i) += h;
            problem->Update(x);
            jacobian.col(i) = (problem->Phi - y0) / h;
        }
        double errJ = (jacobian - J0).norm();
        if (errJ > eps)
        {
            TEST_COUT << "x: " << x0.transpose();
            TEST_COUT << "J*:\n"
                      << jacobian;
            TEST_COUT << "J:\n"
                      << J0;
            TEST_COUT << "(J*-J):\n"
                      << (jacobian - J0);
            ADD_FAILURE() << "Jacobian error out of bounds: " << errJ;
        }
    }
    return true;
}

template <class T>
bool test_jacobian_time_indexed(std::shared_ptr<T> problem, TimeIndexedTask& task, int t, const double eps = 1e-5)
{
    constexpr double h = 1e-6;

    TEST_COUT << "Testing Jacobian with h=" << h << ", eps=" << eps;
    for (int tr = 0; tr < num_trials_; ++tr)
    {
        Eigen::VectorXd x0(problem->N);
        x0 = problem->GetScene()->GetKinematicTree().GetRandomControlledState();
        problem->Update(x0, t);
        TaskSpaceVector y0 = task.Phi[t];
        Eigen::MatrixXd J0 = task.jacobian[t];
        Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(J0.rows(), J0.cols());
        for (int i = 0; i < problem->N; ++i)
        {
            Eigen::VectorXd x = x0;
            x(i) += h;
            problem->Update(x, t);
            jacobian.col(i) = (task.Phi[t] - y0) / h;
        }
        double errJ = (jacobian - J0).norm();
        if (errJ > eps)
        {
            TEST_COUT << "x: " << x0.transpose();
            TEST_COUT << "Phi: " << task.Phi[t].data.transpose();
            TEST_COUT << "J*:\n"
                      << jacobian;
            TEST_COUT << "J:\n"
                      << J0;
            ADD_FAILURE() << "Jacobian error out of bounds: " << errJ;
        }
    }
    return true;
}

UnconstrainedEndPoseProblemPtr setup_problem(Initializer& map, std::string collision_scene = "", std::vector<Initializer> links = std::vector<Initializer>())
{
    Initializer scene;
    if (collision_scene != "")
        scene = Initializer("Scene", {{"Name", std::string("MyScene")}, {"JointGroup", std::string("arm")}, {"Links", links}, {"CollisionScene", std::string(collision_scene)}});
    else
        scene = Initializer("Scene", {{"Name", std::string("MyScene")}, {"JointGroup", std::string("arm")}, {"Links", links}});
    Initializer cost("exotica/Task", {{"Task", std::string("MyTask")}});
    Eigen::VectorXd W = Eigen::Vector3d(3, 2, 1);
    Initializer problem("exotica/UnconstrainedEndPoseProblem", {
                                                                   {"Name", std::string("MyProblem")},
                                                                   {"PlanningScene", scene},
                                                                   {"Maps", std::vector<Initializer>({map})},
                                                                   {"Cost", std::vector<Initializer>({cost})},
                                                                   {"W", W},
                                                               });
    Server::Instance()->GetModel("robot_description", urdf_string_, srdf_string_);

    UnconstrainedEndPoseProblemPtr problem_ptr = std::static_pointer_cast<UnconstrainedEndPoseProblem>(Setup::CreateProblem(problem));

    // Create and test a problem with multiple cost terms
    problem = Initializer("exotica/UnconstrainedEndPoseProblem", {
                                                                     {"Name", std::string("MyProblem")},
                                                                     {"PlanningScene", scene},
                                                                     {"Maps", std::vector<Initializer>({map})},
                                                                     {"Cost", std::vector<Initializer>({cost, cost, cost})},
                                                                     {"W", W},
                                                                 });

    test_random(std::static_pointer_cast<UnconstrainedEndPoseProblem>(Setup::CreateProblem(problem)));

    return problem_ptr;
}

UnconstrainedTimeIndexedProblemPtr setup_time_indexed_problem(Initializer& map)
{
    Initializer scene("Scene", {{"Name", std::string("MyScene")}, {"JointGroup", std::string("arm")}});
    Initializer cost("exotica/Task", {{"Task", std::string("MyTask")}});
    Eigen::VectorXd W = Eigen::Vector3d(3, 2, 1);

    Initializer problem("exotica/UnconstrainedTimeIndexedProblem", {{"Name", std::string("MyProblem")},
                                                                    {"PlanningScene", scene},
                                                                    {"Maps", std::vector<Initializer>({map})},
                                                                    {"Cost", std::vector<Initializer>({cost})},
                                                                    {"W", W},
                                                                    {"T", std::string("10")},
                                                                    {"tau", std::string("0.05")}});
    Server::Instance()->GetModel("robot_description", urdf_string_, srdf_string_);

    UnconstrainedTimeIndexedProblemPtr problem_ptr = std::static_pointer_cast<UnconstrainedTimeIndexedProblem>(Setup::CreateProblem(problem));

    // Create and test a problem with multiple cost terms
    problem = Initializer("exotica/UnconstrainedTimeIndexedProblem", {{"Name", std::string("MyProblem")},
                                                                      {"PlanningScene", scene},
                                                                      {"Maps", std::vector<Initializer>({map})},
                                                                      {"Cost", std::vector<Initializer>({cost, cost, cost})},
                                                                      {"W", W},
                                                                      {"T", std::string("10")},
                                                                      {"tau", std::string("0.05")}});

    test_random(std::static_pointer_cast<UnconstrainedTimeIndexedProblem>(Setup::CreateProblem(problem)));

    return problem_ptr;
}

TEST(ExoticaTaskMaps, testPoint2Line)
{
    try
    {
        TEST_COUT << "PointToLine Test";
        Initializer map("exotica/PointToLine", {{"Name", std::string("MyTask")},
                                                // {"EndPoint", std::string("0.5 0.5 1")},
                                                {"EndPoint", std::string("0.5 0.5 0")},
                                                {"EndEffector", std::vector<Initializer>(
                                                                    {Initializer("Frame", {{"Link", std::string("endeff")},
                                                                                           {"LinkOffset", std::string("0.5 0 0.5")},
                                                                                           {"Base", std::string("base")},
                                                                                           {"BaseOffset", std::string("0.5 0.5 0")}})})}});
        UnconstrainedEndPoseProblemPtr problem = setup_problem(map);
        EXPECT_TRUE(test_random(problem));
        // TODO: Add test_values

        EXPECT_TRUE(test_jacobian(problem));
    }
    catch (...)
    {
        ADD_FAILURE() << "Uncaught exception!";
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    ros::init(argc, argv, "EXOTica_test_maps");
    int ret = RUN_ALL_TESTS();
    Setup::Destroy();
    return ret;
}
