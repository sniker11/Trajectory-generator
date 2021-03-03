#include "umrob/trajectory_generator.h"
#include <fmt/format.h>
#include <numeric>
#include<math.h>
#include <umrob/polynomial.h>

namespace umrob {

TrajectoryGenerator::Segment::Segment( 
    const Eigen::Vector6d& max_velocity,const Eigen::Vector6d& max_acceleration) 
{
    Segment::max_velocity=max_velocity;
    Segment::max_acceleration=max_acceleration;
    

}

TrajectoryGenerator::TrajectoryGenerator(double time_step) : time_step_{time_step} {
   reset();
}

void TrajectoryGenerator::startFrom(const Eigen::Affine3d& pose) {
  reset();
  last_waypoint_ = pose;
  target_pose_ = pose;
}

void TrajectoryGenerator::addWaypoint(const Eigen::Affine3d& pose,const Eigen::Vector6d& max_velocity,const Eigen::Vector6d& max_acceleration){
//state_ = State::ReachingWaypoint;
 target_pose_=pose;
 currentSegment().max_velocity=max_velocity;
 currentSegment().max_acceleration=max_acceleration;
 setSegmentConstraints(last_waypoint_,target_pose_,currentSegment());
 computeSegmentDuration(currentSegment()); 
 computeSegmentParameters(currentSegment()); 

 
}

void TrajectoryGenerator::addWaypoint( const Eigen::Affine3d& pose,double duration) {
    // state_ = State::ReachingWaypoint;
     target_pose_=pose;
     setSegmentConstraints(last_waypoint_,target_pose_,currentSegment());
     currentSegment().duration=duration;
     computeSegmentParameters(currentSegment()); 
}

TrajectoryGenerator::State TrajectoryGenerator::update() {
     if (state_ == State::TrajectoryCompleted) {
        return state_;
    } else {
        state_ = State::ReachingWaypoint;
    }

    auto evaluate = [this]() {
         // TODO evaluate the polynomial and its derivatives into next_pose_vec,
         // target_velocity_ and target_acceleration_

        Eigen::Vector6d next_pose_vec;
        for (int k = 0; k < 6; k++)
       {    auto& poly = currentSegment().polynomials[k];
            next_pose_vec(k) = poly.evaluate(currentSegment().current_time);
            target_velocity_(k) =poly.evaluateFirstDerivative(currentSegment().current_time);
            target_acceleration_(k)=poly.evaluateSecondDerivative(currentSegment().current_time);
       }
        target_pose_.translation() = next_pose_vec.head<3>();
        target_pose_.linear() =Eigen::Quaterniond::fromAngles(next_pose_vec.tail<3>()).toRotationMatrix();

    };

    evaluate();

    currentSegment().current_time += time_step_;
    if (currentSegment().current_time > currentSegment().duration) {
        // TODO switch to next segment and re-evaluate if needed// si on est dans le dernier segment
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
        current_segment_idx_++;
        state_ = State::WaypointReached;
        last_waypoint_=target_pose_;
       if (current_segment_idx_==segments_.size()) 
        {evaluate();
        state_ = State::TrajectoryCompleted;
        }
    }

    return state_;
   
}

const Eigen::Affine3d& TrajectoryGenerator::targetPose() const {
    return target_pose_;
}

const Eigen::Vector6d& TrajectoryGenerator::targetVelocity() const {
    return target_velocity_;
}

const Eigen::Vector6d& TrajectoryGenerator::targetAcceleration() const {
    return target_acceleration_;
}

TrajectoryGenerator::State TrajectoryGenerator::state() const {
    return state_;
}

double TrajectoryGenerator::duration() const {
 double Duree;
 size_t i;
 for (i=0;i<segments_.size();i++)
 Duree+=segments_.at(i).duration;

return Duree;
}
void TrajectoryGenerator::reset(){
    state_ = State::Idle;
    last_waypoint_.setIdentity();
    segments_.clear();
    current_segment_idx_ = 0;
    target_pose_ = last_waypoint_;
    target_velocity_.setZero();
    target_acceleration_.setZero();

}

void TrajectoryGenerator::setSegmentConstraints(const Eigen::Affine3d& from, const Eigen::Affine3d& to,Segment& segment) {
    Eigen::Vector6d from_vec, to_vec;
    from_vec << from.translation(),Eigen::Quaterniond(from.linear()).getAngles();
    to_vec << to.translation(), Eigen::Quaterniond(to.linear()).getAngles();

    for (size_t i = 0; i < 6; i++) {
        auto& poly = segment.polynomials[i];
        auto& cstr = poly.constraints();
        cstr.xi = 0;
        cstr.xf = 0;
        cstr.yi = from_vec(i);
        cstr.yf = to_vec(i);
    }
}
//! Tmin_vel = 30Δy/16vmax
//! Tmin_acc = sqrt(10sqrt(3)Δy/3amax);
void TrajectoryGenerator::computeSegmentDuration (Segment& segment) {

Eigen::Vector6d Tmin_vel;
Eigen::Vector6d Tmin_acc;
double DP;
double w{0};
Eigen::Vector6d SEGDUR;
    for (int i = 0; i < 6; i++) {
DP=(segment.polynomials[i].constraints().yf)-(segment.polynomials[i].constraints().yi);
Tmin_vel(i)=(30*DP)/(16*segment.max_velocity(i));
Tmin_acc(i)= sqrt((10*sqrt(3)*(DP)/(3*segment.max_acceleration(i))));
SEGDUR(i)=std::max(Tmin_acc(i),Tmin_vel(i));
if (SEGDUR(i)>w) w=SEGDUR(i);
    }
 segment.duration=w; 
}
 

void TrajectoryGenerator::computeSegmentParameters(Segment& segment) {

    for (size_t i=0;i<6;i++)
    {
        auto& poly = segment.polynomials[i];
        poly.constraints().xf=segment.duration;
        poly.computeCoefficients();
    }

    //to do implement
}

TrajectoryGenerator::Segment& TrajectoryGenerator::currentSegment(){
    return segments_.at(current_segment_idx_);
}

} // namespace umrob
