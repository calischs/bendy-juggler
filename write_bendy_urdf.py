#!/usr/bin/env python
from __future__ import division
import sys
import xml.etree.cElementTree as ET

def indent(elem, level=0):
    '''
    http://stackoverflow.com/questions/749796/pretty-printing-xml-in-python
    '''
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def write_bendy(l,m,extra_mass,k,c,n_links,filename='bendy.xml'):
    l = float(l); m = float(m); extra_mass = float(extra_mass)
    k = float(k); c = float(c); n_links = int(n_links);
    robot = ET.Element("robot")
    robot.set("xmlns","http://drake.mit.edu")
    robot.set("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance")
    robot.set("xsi:schemaLocation","http://drake.mit.edu ../../doc/drakeURDF.xsd")
    robot.set("name","bendy")

    baselink = ET.SubElement(robot, "link")
    baselink.set("name", "baselink")
    viz = ET.SubElement(baselink,"visual")
    geo = ET.SubElement(viz,"geometry")
    box = ET.SubElement(geo,"box")
    box.set("size",".05 .05 .05")
    mat = ET.SubElement(viz,"material")
    mat.set("name","green")
    col = ET.SubElement(mat,"color")
    col.set("rgba","0 1 0 1")

    link_len = l/n_links
    print link_len
    link_mass = m/n_links
    for i in range(n_links):
        link = ET.SubElement(robot, "link")
        link.set("name", "link%02d"%(i+1))
        
        inertial = ET.SubElement(link,"inertial")
        orig = ET.SubElement(inertial,"origin")
        orig.set("xyz","0 0 %.3f"%(-(i)*link_len))
        orig.set("rpy","0 0 0")
        mass = ET.SubElement(inertial,"mass")
        mass.set("value","%.3f"%link_mass)
        inertia = ET.SubElement(inertial,"inertia")
        inertia.set("ixx","1.")
        inertia.set("ixy","0.")
        inertia.set("ixz","0.")
        inertia.set("iyy",".33")
        inertia.set("iyz","0.")
        inertia.set("izz","1.")
        
        viz = ET.SubElement(link,"visual")
        orig = ET.SubElement(viz,"origin")
        orig.set("xyz","0 0 %.3f"%(-.5*link_len)) #set geometry in link coords
        orig.set("rpy","0 0 0")
        geo = ET.SubElement(viz,"geometry")
        cyl = ET.SubElement(geo,"cylinder")
        cyl.set("length","%.3f"%(link_len))
        cyl.set("radius","%.3f"%(.01))
        mat = ET.SubElement(viz,"material")
        mat.set("name","link%02dcol"%(i+1))
        col = ET.SubElement(mat,"color")
        col.set("rgba","0 %.2f %.2f 1"%(i/(n_links-1.),1.-i/(n_links-1.)))
        
        '''
        #we got rid of collisions for the hybrid model
        coll = ET.SubElement(link,"collision")
        geo = ET.SubElement(coll,"geometry")
        cyl = ET.SubElement(geo,"cylinder")
        cyl.set("length","%.3f"%(link_len))
        cyl.set("radius","%.3f"%(.05))
        '''

    for i in range(n_links):
        joint = ET.SubElement(robot, "joint")
        joint.set("name", "joint%02d"%(i+1))
        joint.set("type", "continuous")
        parent = ET.SubElement(joint, "parent")
        parent.set("link","baselink" if i==0 else "link%02d"%(i))
        child = ET.SubElement(joint, "child")
        child.set("link","link%02d"%(i+1))
        orig = ET.SubElement(joint, "origin")
        orig.set("xyz","0 0 %.3f"%(0 if i==0 else -link_len))
        axis = ET.SubElement(joint, "axis")
        axis.set("xyz","0 1 0")
        dynam = ET.SubElement(joint, "dynamics")
        dynam.set("damping","%.3f"%(c))
    for i in range(1,n_links):
        force_elem = ET.SubElement(robot, "force_element")
        force_elem.set("name", "joint_spring%d"%(i))
        spring = ET.SubElement(force_elem, "torsional_spring")
        spring.set("stiffness","%.3f"%k)
        spring.set("rest_angle","0.")
        joint_ref = ET.SubElement(spring, "joint")
        joint_ref.set("name","joint%02d"%(i+1))
    
    ### HACK
    #force_elem = ET.SubElement(robot, "force_element")
    #force_elem.set("name", "joint_spring%d"%(0))
    #spring = ET.SubElement(force_elem, "torsional_spring")
    #spring.set("stiffness","%.3f"%1e4)
    #spring.set("rest_angle","-1.57")
    #joint_ref = ET.SubElement(spring, "joint")
    #joint_ref.set("name","joint%02d"%(0+1))
    ### HACK

    link = ET.SubElement(robot, "link")
    link.set("name", "final_mass")
    viz = ET.SubElement(link,"visual")
    orig = ET.SubElement(viz,"origin")
    orig.set("xyz","0 0 %.3f"%(-link_len)) #set geometry in link coords
    geo = ET.SubElement(viz,"geometry")
    sphere = ET.SubElement(geo,"sphere")
    sphere.set("radius","%.3f"%(.02))
    mat = ET.SubElement(viz,"material")
    mat.set("name","mass")
    col = ET.SubElement(mat,"color")
    col.set("rgba","0 0 1 1")
    inertial = ET.SubElement(link,"inertial")
    orig = ET.SubElement(inertial,"origin")
    orig.set("xyz","0 0 %.3f"%(-link_len))
    mass = ET.SubElement(inertial,"mass")
    mass.set("value","%.3f"%extra_mass)
    joint = ET.SubElement(robot, "joint")
    joint.set("name", "mass_weld")
    joint.set("type","fixed")
    parent = ET.SubElement(joint, "parent")
    parent.set("link","link%02d"%(n_links))
    child = ET.SubElement(joint, "child")
    child.set("link","final_mass")
    
    trans = ET.SubElement(robot, "transmission")
    trans.set("type","SimpleTransmission")
    trans.set("name","base_trans")
    act = ET.SubElement(trans, "actuator")
    act.set("name","base_torque")
    joint = ET.SubElement(trans, "joint")
    joint.set("name","joint01")
    red = ET.SubElement(trans, "mechanicalReduction")
    red.text = "1"
    
    frame = ET.SubElement(robot, "frame")
    frame.set("name","end_effector")
    frame.set("link","link%02d"%(n_links))
    frame.set("xyz","0 0 %.3f"%(-link_len))
    frame.set("rpy","0 0 0")


    tree = ET.ElementTree(robot)
    indent(robot)
    tree.write(filename)


if __name__ == '__main__':
	#write_bendy(l,m,extra_mass,k,n_links,filename='bendy.xml'):
	write_bendy(*sys.argv[1:])



