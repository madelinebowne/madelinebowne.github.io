#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:39:05 2024

@author: madelinebowne
"""

class WarehouseOld:
    def __init__(self, env, id, altitude, e, i, RAAN_0, omega, theta_0, plane_location, dry_mass, payload_mass, fuel_mass, ISP, available_slot, spare_supply, status, initial_capacity, max_capacity,satellites,planes,numSpares,drift_fuel_per_second,fuel_cost_per_kg,adr_operation_cost,altitude_customer,altitude_parking,sat_dry_mass,sat_fuel_mass,sat_ISP,numSatellitesPerPlane,numPlanes, planeSpot, ADRassignment,customer,servicer,launcher,fuel_mass_wh,satellite_drymass,satellite_fuel,electric_fuel_price_per_kg,chemical_fuel_price_per_kg,totalTime,numADR,ADR_dry_mass,ADR_fuel_mass):
        self.env = env
        self.id = id
        self.altitude = altitude  # km
        self.e = e  # Assuming circular orbits
        self.i = i  # Inclination in degrees
        self.RAAN_0 = RAAN_0  # Initial RAAN
        self.omega = omega  # Argument of periapsis for circular orbits
        self.theta_0 = theta_0  # Initial true anomaly
        self.plane_location = plane_location  # Dictates which plane index it's in, 0 if not
        self.dry_mass = dry_mass  # Dry mass of warehouse alone
        self.payload_mass = payload_mass  # Mass of satellites onboard
        self.fuel_mass = fuel_mass  # How much mass in fuel
        self.ISP = ISP
        self.spare_supply = spare_supply  # How many new satellites aboard
        self.available_slot = available_slot  # How many available slots for old satellites
        self.status = status  # 'available' or 'unavailable'
        self.max_capacity = max_capacity
        self.initial_capacity = initial_capacity
        self.container = simpy.Container(env, capacity=max_capacity, init=initial_capacity)
        self.collectedSatellites = 0  # List of satellites stored in the warehouse
        self.satellites = []
        self.propagate_time_tracker = 0  # Epoch time (starting time)
        self.NoSpares = False
        self.resupplyOnTheWay = False
        self.criticalSpares = False
        self.criticalSpareLevel = round(0.61*self.max_capacity)
        self.noSpare_event = env.event()
        self.below_critical_spare_event = env.event()
        self.resupply_event = env.event()
        self.low_fuel_event = env.event()
        self.adr_level_event = env.event()
        self.sat_delivery_event = env.event()
        self.costs = 0
        self.upgrade2refuelsats = 0 #boolean is 1
        self.upgrade2repairsats = 0 #boolean is 1
        self.refuelUpgraded = 0 #becomes 1 once wh is able to refuel collected sats
        self.repairUpgraded = 0 #becomes 1 once wh is able to repair collected sats
        self.upgrade2repairCost = 0 #based on uncertain function
        self.upgrade2refuelCost = 0 #based on uncertain function
        self.upgrade2refuelTime = 1000 #amount of time (on top of typical launch time) to include refuel module in payload
        self.upgrade2repairTime = 1000 #amount of time (on top of typical launch time) to include repair module in payload
        self.deorbit = 0 #if servicer decides to de-orbit the entire warehouse, this will be set to 1
        self.planeTracker = planeSpot #the last plane the warehouse drifted to. Starts at 0 for initialization
        self.ADRtracker = ADRassignment #contains the array of ID of the attached ADR vehicles, 0 if none
        #self.ADRtracker = self.create_array(ADRassignment, numADR)
        self.adrtargetlevel = 0
        self.numResupply = 0
        self.wh_fuel_for_sats = 100 #total fuel depot for warehouse to refuel satellites
        self.fuel_for_sats = 0 #how much fuel is left to refuel satellites

        #this schedules warehouse resupply when it reaches its critical container level, accounting for satellites it has collected
        self.env.process(self.warehouseResupply(env,customer,servicer,launcher,satellites,fuel_mass_wh,satellite_drymass,satellite_fuel,electric_fuel_price_per_kg,chemical_fuel_price_per_kg,drift_fuel_per_second,adr_operation_cost,altitude_customer,altitude_parking,sat_dry_mass,sat_ISP,numSatellitesPerPlane,ADR_dry_mass,ADR_fuel_mass))

        #this process has the warehouse drift between orbital planes and release spare/active satellites as needed when it reaches each plane  
        self.env.process(self.drift_and_replenish(env,satellites,planes,numSpares,drift_fuel_per_second,fuel_cost_per_kg,adr_operation_cost,altitude_customer,altitude_parking,sat_dry_mass,sat_fuel_mass,sat_ISP,numSatellitesPerPlane,numPlanes,totalTime,customer,servicer))
    
    def create_array(self, ADRassignment, numADR):
        # Create an array of zeros with length Y
        arr = np.zeros(numADR, dtype=int)
        
        # Set the value at index X to 1
        if ADRassignment <= numADR:
            arr[(ADRassignment-1)] = 1
        else:
            raise ValueError(f"ADRassignment {ADRassignment} should be less than numADR {numADR}")
        
        return arr
        
    def propagate(self, env):
        deltaT = env.now - self.propagate_time_tracker
        SEA_LEVEL_GRAV_ACCEL = 9.81
        EARTH_RADIUS = pk.EARTH_RADIUS/1000 #in km
        J2 = pk.EARTH_J2
        mu = (pk.MU_EARTH/1000000000)
        #altitude must be in semi-major axis form
        
        if deltaT > 0:
            v = math.sqrt(SEA_LEVEL_GRAV_ACCEL * EARTH_RADIUS ** 2 / ((self.altitude)+EARTH_RADIUS))  # Calculate the orbital velocity [m/s] based on the current radius
            n = v / ((self.altitude)+EARTH_RADIUS)  # Calculate the mean motion [rad/second], which is the rate of orbit
            # Calculate J2 perturbations
            delta_drift = compute_delta_drift(((self.altitude)+EARTH_RADIUS) , self.i)
            # Update RAAN
            self.RAAN_0 = (self.RAAN_0 + np.rad2deg(abs(delta_drift * deltaT))) % 360

            # Calculate positions (assuming simple propagation for demonstration)
            theta_dot = np.sqrt(mu / (((self.altitude)+EARTH_RADIUS)  ** 3)) * deltaT  # Mean motion * time
            self.theta_0 = (self.theta_0 + abs(theta_dot)) % 360

            # Update last propagated time
            self.propagate_time_tracker = env.now # Update the epoch time to the current simulation time
            
            #print(f"warehouse {self.id} is now at RAAN {self.RAAN_0} after drifting {abs(delta_drift * deltaT)} degrees")
    

    
    def release_new_sat(self,env,orbitalSlot,closest_plane,RAAN,theta,status,satellites,drift_fuel_per_second,fuel_cost_per_kg,adr_operation_cost,altitude_customer,altitude_parking,sat_dry_mass,sat_fuel_mass,sat_ISP,numSatellitesPerPlane,numPlanes,totalTime,customer,servicer,planes):
        """
        The warehouse container level includes both new spares and collected satellites. To release a new satellite, make sure that the container level contains more than just old sats.
        """
        if (self.container.level-self.collectedSatellites) > 0:
            
            #print(f"inside release_new_sat, the RAAN is {RAAN} and self RAAN is {self.RAAN_0}")
            """update the location of both the warehouse and the closest_plane, which is provided as an input when the warehouse has drifted to it"""
            self.propagate(env)
            closest_plane.propagate(env)
            
            #print(f"release_new_sat now turned on for wh {self.id} with RAAN {self.RAAN_0} and plane {closest_plane.id} with RAAN {closest_plane.RAAN_0}")
            
            """orbital_maneuver_special is an adjusted version of orbital_maneuver and it only considers in-plane maneuvering.
            This version sets self equal to the chaser, even though the warehouse is not actually moving to the failed/retired sat, but releasing a sat that will maneuver to that location.
            Not that CBmode=1, which means that the mass and ISP for the chaser is manually input rather than coming from the warehouse object properties, since the warehouse is not actually moving.
            The warehouse is only listed as the chaser to extract its position information.
            This assumes that the released satellite has high-impulse thrust (chemical fuel)"""
            
            #print(f"the self raan is {self.RAAN_0} and self theta is {self.theta_0} and self altitude is {self.altitude} and RAAN is {RAAN}, the theta is {theta} and {altitude_customer} and {altitude_parking}")
            
            rend_cost_wh, total_fuel_wh, total_time_wh, deltaV_wh, man_mode_wh = orbital_maneuver_special(env,self,RAAN,theta,altitude_customer,altitude_parking,fuel_cost_per_kg,drift_fuel_per_second,1,sat_dry_mass,sat_ISP)
            yield env.timeout(total_time_wh)
            
            #check what current active satellite situation is and make sure there aren't more than there should be
            #actives_in_plane = [active for active in satellites if active.plane_location == closest_plane.id and active.status == 'active']
            #if len(actives_in_plane) >= numSatellitesPerPlane:
            #    status = 'spare'
            
            
            """when a satellite leaves the warehouse, the warehouse payload mass is reduced by its total mass (dry+fuel)"""
            self.payload_mass -= (sat_dry_mass+sat_fuel_mass)
            yield self.container.get(1)
            
            
            #print(f"warehouse {self.id} will release new sat for spare {sat.id} in plane {sat.plane_location} in time {self.env.now+total_time_wh}")
            
            # Simulate the time it takes for the spare satellite to drift to the required location
            
            """wait for the time it will take the new satellite to reach the position of the failed/retired satellite, then propagate"""
            #yield env.timeout(total_time_wh)
            self.propagate(env)
            for sat in satellites:
                sat.propagate(env)
            closest_plane.propagate(env)
            closest_plane.incomingSpares = 0
            closest_plane.incomingActives = 0
            
            #check what current active satellite situation is and make sure there aren't more than there should be
            actives_in_plane = [active for active in satellites if active.plane_location == closest_plane.id and active.status == 'active']
            if len(actives_in_plane) >= numSatellitesPerPlane:
                status = 'spare'
                #warehouse delivery updates
            
            #get the RAAN and theta from the updated plane properties
            """get the RAAN from the plane properties and adjust the given theta for the time passed (based on orbital slot in need from drift_and_replenish)"""
            RAAN_0 = closest_plane.RAAN_0
            theta_dot = np.sqrt(((pk.MU_EARTH/1000000000)) / (((altitude_customer)+(pk.EARTH_RADIUS/1000))  ** 3)) * total_time_wh  # Mean motion * time
            theta_0 = (theta + abs(theta_dot)) % 360
            
            #theta_0 = closest_plane.orbitalSlots[(orbitalSlot-1), 1]
            
            """"get a new failNum for the new satellite. Adjusted to consider the remaining amount of time left. Added to present time to have present-time failure events"""
            newTime = totalTime-env.now
            #print(f"the new time is {newTime}")
            failArray = assign_failure(1,int(newTime))+env.now
            
            """the new id for the sat is one higher than the highest existing id"""
            max_id = max(sat.id for sat in satellites)
            
            satellite_properties = satellites[max_id-1]
    
            """Create a new satellite object with a new ID and univeral properties extracted from the most recent satellite
            RAAN and theta are from updated values"""
                
            
            new_satellite_id = max_id + 1    
            new_lifeLeft = satellite_properties.design_life
            new_satellite = Satellite(
                env,
                new_satellite_id,
                satellite_properties.altitude,
                satellite_properties.e,
                satellite_properties.i,
                RAAN_0,
                0,  # Assuming initial mean anomaly is 0 for simplicity
                theta_0,
                closest_plane.id,
                0,  # warehouse_location
                sat_dry_mass,
                sat_fuel_mass,
                status,
                new_lifeLeft,
                satellite_properties.servicePrice,
                satellite_properties.design_life,
                satellite_properties.failurePenalty,
                satellite_properties.failureRate,
                failArray[0],
                orbitalSlot,
                customer,servicer,totalTime,planes
            )
            satellites.append(new_satellite)
            #print(f"Satellite in plane {orbitalSlot} with RAAN {closest_plane.RAAN_0} replaced by sat {new_satellite.id} with status {new_satellite.status} from warehouse {self.id} with RAAN {self.RAAN_0} at time {self.env.now}")
            #sat_fuel_mass * (new_lifeLeft / satellite_properties.design_life) 
        
        
                            # Trigger status event for the new satellite
            if status == 'spare':
                closest_plane.spareNum += 1
                new_satellite.activation_event = env.event()
            if status == 'active':
                closest_plane.activeNum += 1
                if not closest_plane.activeNumChange_event.triggered:
                    print("activeNumChange has been triggered")
                    closest_plane.activeNumChange_event.succeed()
                new_satellite.activation_event = env.event()
                new_satellite.activation_event.succeed()
            if not self.sat_delivery_event.triggered:
                self.sat_delivery_event.succeed()
            #self.sat_delivery_event = env.event()
            new_satellite.propagate_time_tracker = self.env.now  # Epoch time (starting time)
            new_satellite.retire_time_tracker = self.env.now
            new_satellite.failed = False
            new_satellite.retired = False
            new_satellite.old = False
            new_satellite.failure_event = env.event()
            new_satellite.retirement_event = env.event()
            new_satellite.old_event = env.event()
            new_satellite.criticalValue = satellite_properties.criticalValue #2.628e+6 #one month in seconds
            #new_satellite.failNum = 1*10**100
            new_satellite.orbitalSlotIndex = orbitalSlot   
            #new_satellite.scheduledRetirement = satellite_properties.design_life
            #new_satellite.scheduledNearRetirement = satellite_properties.design_life-satellite_properties.criticalRetirementTime
            #new_satellite.scheduledDeorbit = satellite_properties.design_life+satellite_properties.criticalValue

            
            """if by releasing this satellite, the warehouse reaches a critical number of spares (parametric variable), it triggers an event that will request a resupply to the warehouse"""
            if (self.container.level-self.collectedSatellites) <= self.criticalSpareLevel:
                self.criticalSpares = True
                if not self.below_critical_spare_event.triggered:
                    #print(f"warehouse {self.id} has reached the critical spare level")
                    self.below_critical_spare_event.succeed()
                    
            if (self.max_capacity-self.container.level) <= self.adrtargetlevel:
                if not self.adr_level_event.triggered:
                    self.adr_level_event.succeed()
                    self.adr_level_event = env.event()
                    self.adrtargetlevel = 0
                    
                
        #if there are no spares available, it triggers a no spare event
        else:
            #print(f"warehouse couldnt complete replacement request for sat {sat.id} because its container is empty")
            self.NoSpares = True
            if not self.noSpare_event.triggered:
                self.noSpare_event.succeed()
        
    
    def drift_and_replenish(self,env,satellites,planes,numSpares,drift_fuel_per_second,fuel_cost_per_kg,adr_operation_cost,altitude_customer,altitude_parking,sat_dry_mass,sat_fuel_mass,sat_ISP,numSatellitesPerPlane,numPlanes,totalTime,customer,servicer):
        while True:
            #print(f"spare event is {self.below_critical_spare_event}")
            """drift_and_replenish moves the warehouse through each orbital spare. As it crosses each plane, it looks to see how many spares and actives need to be added to the plane.
            It calls release_new_sat for each replacement spare, starting by propagating itself, each plane, and each satellite"""
            self.propagate(env)
            for plane in planes:
                plane.propagate(env)
            for sat in satellites:
                sat.propagate(env)
            
            #print(f"the current wh plane for wh {self.id} is {self.planeTracker}")
            
            """it figures out which plane is the next to drift to. In the event a warehouse is initialized between planes, it rounds to the nearest plane.
            If it starts directly in an integer plane, it simply moves to the next integer plane."""
            if self.planeTracker % 1 != 0:
                #print(f"the plane track is {(math.ceil(self.planeTracker))}")
                upNext = (math.ceil(self.planeTracker))
                if upNext != numPlanes:
                    nextPlane = (math.ceil(self.planeTracker)) % (numPlanes)
                else:
                    nextPlane = upNext
                #print(f"the next plane is {(math.ceil(self.planeTracker)) % (numPlanes)}")
            elif self.planeTracker % 1 == 0:
                upNext = (self.planeTracker+1)
                if upNext != numPlanes:
                    nextPlane = (self.planeTracker+1) % (numPlanes)
                else:
                    nextPlane = upNext

            #get the object associated with the particular plane index
            plane_closest_list = [plane for plane in planes if plane.id == nextPlane] 
            plane_closest = plane_closest_list[0]
            
            #print(f"the next wh plane for wh {self.id} is {plane_closest.id}")
            
            """calculate the amount of time it will take for the plane and warehouse to drift to the same RAAN, taking their altitude and initial RAAN into account"""
            #print(f"warehouse {self.id} at RAAN {self.RAAN_0} and the plane_closest is at {plane_closest.RAAN_0}")
            DriftTime = time_to_drift_same_plane(self.RAAN_0, plane_closest.RAAN_0, self.altitude, altitude_customer, self.i, self.i)
            
            #print(f"the next wh plane for wh {self.id} is {plane_closest.id} and it will take {DriftTime} seconds with self.RAAN_0 {self.RAAN_0} and plane_closest.RAAN_0 {plane_closest.RAAN_0}")
            
            #print(f"the next wh plane for wh {self.id} is {plane_closest.id} in time {DriftTime}")
            
            #print(f"before yeild, in drift and replen, self.RAAN_0 is {self.RAAN_0} and plane RAAN is {plane_closest.RAAN_0}")
            
            """propagate the warehouse, planes, and satellites for that designated amount of time"""
            yield self.env.timeout(DriftTime)
            self.propagate(env)
            #plane_closest.propagate(env)
            for plane in planes:
                plane.propagate(env)
            for sat in satellites:
                sat.propagate(env)
            
            #print(f"the current operational wh plane (RAAN {self.RAAN_0} for wh {self.id} is {plane_closest.id} with RAAN {plane_closest.RAAN_0}")
            
            #print(f"after drift time {DriftTime}, warehouse {self.id} has RAAN {self.RAAN_0} is now in plane {plane_closest.id} that has RAAN {plane_closest.RAAN_0}")
            #update the warehouse planeTracker to reflect its current plane location 
            self.planeTracker = plane_closest.id

            
            #print(f"in drift and replen, self.RAAN_0 is {self.RAAN_0} and plane RAAN is {plane_closest.RAAN_0}")
            
            #figure out how many active and spare satellites are in that plane
            activeSatsss = [sat for sat in satellites if sat.status == 'active' and sat.plane_location == plane_closest.id] #not assigned to be collected, not replaced
            spareSatsss = [sat for sat in satellites if sat.status == 'spare' and sat.plane_location == plane_closest.id]
            
            statutus = [sat.status for sat in satellites if sat.plane_location == plane_closest.id]
            
            #activeIds = [sat.plane_location for sat in satellites if sat.status == 'active' and sat.assignment == 'unassigned']
            
            #figure out how many actives and spares the warehouse needs to release based on how mant should be there
            #activesInNeed = numSatellitesPerPlane - len(activeSatsss)
            
            #make sure another warehouse didn't already initiate a resupply to the orbital plane
            if plane_closest.incomingSpares == 0:
                activesInNeed = numSatellitesPerPlane - len(activeSatsss)
                sparesInNeed = numSpares - len(spareSatsss)
                plane_closest.incomingActives = activesInNeed
                plane_closest.incomingSpares = sparesInNeed
                #if activesInNeed > 0 or sparesInNeed > 0:
                #    plane_closest.incomingSpares = 1
            else:
                activesInNeed = 0
                sparesInNeed = 0
                
            #sparesInNeed = numSpares - len(spareSatsss)
            
            #print(f"actives in need are {activesInNeed}")
            #print(f"spares in need are {sparesInNeed}")
            #print(f"the satellites active are {activeIds}")
            
            #print(f"the empty orbital slots in plane {plane_closest.id} are {plane_closest.orbitalSlotOccupancy}")
            #print(f"the spares in need are {sparesInNeed} and the number of spare sats is {len(spareSatsss)} long")
            #print(f"the actives in need are {activesInNeed} and the number of actives sats is {len(activeSatsss)} long")
            
            RAAN = plane_closest.RAAN_0
            #print(f"inside dirft and replenish, RAAN is {RAAN}")

            #print(f"in plane {self.planeTracker} there are {activesInNeed} actives in need and {sparesInNeed} spares in need")
            
            """to assign the satellites to the correct theta values (orbitalSlots) look at the plane's orbitalSlot variable.
            Figure out which slots are designated for active satellites and which are for spare satellites.
            Figure out how many of those slots are occupied, openSlots gives orbitalSlots that don't have spare and active satellites.
            openSlots_excludeSpares give just the empty slots designated for active sats"""
            totalSlots = np.arange(1, len(plane_closest.orbitalSlots[:,1])+1)
            #spareSlots = np.arange(numSpares+1, len(plane_closest.orbitalSlots[:,1])+1)
            spareSlots = np.arange(numSatellitesPerPlane+1, len(plane_closest.orbitalSlots[:,1])+1)
            occupiedSlots = [sat.orbitalSlotIndex for sat in satellites if (sat.status == 'active' or sat.status == 'spare') and sat.plane_location == plane_closest.id]
            #occupiedSlots = [sat.orbitalSlotIndex for sat in satellites if (sat.status == 'active' or sat.status == 'spare') and sat.plane_location == plane_closest.id and sat.assignment == 'unassigned' and sat.replacementStatus == 0]
            
            occupiedSlots_active = [sat.orbitalSlotIndex for sat in satellites if sat.status == 'active' and sat.plane_location == plane_closest.id]
            occupiedSlots_spare = [sat.orbitalSlotIndex for sat in satellites if sat.status == 'spare' and sat.plane_location == plane_closest.id]
            
            
            #print(f"occupiedSlots_spare are {occupiedSlots_spare} and occupiedSlots_active are {occupiedSlots_active}")
           
            openSlots = np.setxor1d(occupiedSlots, totalSlots)
            #openSlots_excludeSpares = np.setxor1d(openSlots, spareSlots)
            excludeSpares = np.setxor1d(totalSlots, occupiedSlots_spare)
            excludeActives= np.setxor1d(totalSlots, occupiedSlots_active)
            openSlots2 = np.setxor1d(excludeSpares,occupiedSlots_active)
            openSlots3 = np.setxor1d(excludeActives,occupiedSlots_spare)
            
            openSlots_excludeSpares  = openSlots2
            
            #print(f" the totalSlots is {totalSlots} with len {len(totalSlots)} and spareSlots is len{len(spareSlots)} and occupiedSlots is {occupiedSlots} with len {len(occupiedSlots)} and openSlots is {openSlots} with len {len(openSlots)}")
            #print(f"the status for all sats in this plane {plane_closest.id} are {statutus}, and len(plane_closest.orbitalSlots[:,1])+1 is {len(plane_closest.orbitalSlots[:,1])+1} and numSatellitesPerPlane is {numSatellitesPerPlane} and len(activeSatsss) is {len(activeSatsss)} and openSlots_excludeSpares is {len(openSlots_excludeSpares)} and range(activesInNeed) is {range(activesInNeed)}")
            
            #6 spare, 11 active, 46 total slots, 31 open slots
            
            #print(f"totalSlots are {totalSlots}, spareSlots are {spareSlots}, occupiedSlots are {occupiedSlots}, openSlots are {openSlots}, openSlots_excludeSpares are {openSlots_excludeSpares} ")
            
            #print(f"the openSlots_excludeSpares is {openSlots_excludeSpares} and activesInNeed are {activesInNeed}")
            
            #start with replenishing actives, if necessary
            #orbitalSlotActive = 1
            
            """if actives need to be delivered, it loops through each need.
            Index is the slot number associated with each empty active slot.
            It uses this index to get the associate theta value.
            the orbitalSlotActive number is simply the "slot number" associated with the slot, which will be saved in the satellite's properties.
            The status is set to immediately go to active, and then release_new_sat is called with these specific properties"""
            
            if activesInNeed > 0:
                for slot in range(activesInNeed):
                ##for slot in range(len(openSlots_excludeSpares)):
                    #print(f"the slot is {slot} and activesInNeed is {activesInNeed} and range(activesInNeed) is {range(activesInNeed)} and openSlots_excludeSpares is {openSlots_excludeSpares} and openSlots_excludeSpares length is {len(openSlots_excludeSpares)}")
                    #print(f"the slot is {slot}")
                    index = int(openSlots_excludeSpares[slot])-1 #this -1 was added to fix a bug but I'm not sure this is actually correct
                    theta = plane_closest.orbitalSlots[index,1]
                    orbitalSlotActive  = index+1
                    #orbitalSlotActive = 0
                    status = 'active'
                    orbitalSlot = orbitalSlotActive
                    ##print(f"plane {plane_closest.id} at RAAN {plane_closest.RAAN_0} is missing a spare sat, so wh {self.id} in plane {self.planeTracker} with RAAN {self.RAAN_0} is releasing a new sat to slot {index}")
                    env.process(self.release_new_sat(env,orbitalSlot,plane_closest,RAAN,theta,status,satellites,drift_fuel_per_second,fuel_cost_per_kg,adr_operation_cost,altitude_customer,altitude_parking,sat_dry_mass,sat_fuel_mass,sat_ISP,numSatellitesPerPlane,numPlanes,totalTime,customer,servicer,planes))
                    ##orbitalSlotActive += 1
                    
            """is spares need to be delivered, it loops through each need.
            Spares are also activated in order, so no need to get an index -- it knows which slots are empty based on how many spares need to be delivered.
            the status is set to spare and then release_new_sat is called with its specific properties"""
            #now replenish spares
            orbitalSlotSpare = 1
            if sparesInNeed > 0:
                for slot in range(sparesInNeed):
                    
                    #print(f"the index is {index} and number of sparesInNeed is {len(sparesInNeed)}")
                    #theta = plane_closest.orbitalSlots[index-1, 1]
                    #theta = plane_closest.orbitalSlots[slot+numSpares,1]
                    #theta = 1
                    index = int(openSlots3[slot])-1
                    theta = plane_closest.orbitalSlots[index,1]
                    orbitalSlotSpare  = index+1
                    orbitalSlot = orbitalSlotSpare
                    #theta = plane_closest.orbitalSlots[slot+numSpares,1]
                    status = 'spare'
                    #orbitalSlot = orbitalSlotSpare
                    #print(f"plane {plane_closest.id} at RAAN {plane_closest.RAAN_0} is missing a spare sat, so wh {self.id} in plane {self.planeTracker} with RAAN {self.RAAN_0} is releasing a new sat to slot {index}")
                    env.process(self.release_new_sat(env,orbitalSlot,plane_closest,RAAN,theta,status,satellites,drift_fuel_per_second,fuel_cost_per_kg,adr_operation_cost,altitude_customer,altitude_parking,sat_dry_mass,sat_fuel_mass,sat_ISP,numSatellitesPerPlane,numPlanes,totalTime,customer,servicer,planes))
                    orbitalSlotSpare += 1
    

    def receive_satellite(self, env, adr,num_satellites, satellites,IDs,target_level):
        """if a warehouse receives an old satellite from an ADR delivery, this triggers for the number of satellite in the collection.
        First, it checks that it has spare to receive the satellites. If it does, it uses container.put, adds the satellite IDs to its object properties,
        appends the satellite to its list of objects, and updates it payload mass to include the sat dry mass for each satellite
        If there isn't enough room, monitor_collection will wait until the warehouse releases enough satellites to make room"""
        #print(f"Warehouse {self.id} is getting {num_satellites} satellites and has container level of {self.container.level} at time {self.env.now}.")
        if (self.max_capacity-self.container.level) >= num_satellites:
            print(f"Warehouse {self.id} to receive {num_satellites} satellite(s) at time {self.env.now}")
            yield self.container.put(num_satellites)
            self.satellites.append(IDs)
            self.collectedSatellites += num_satellites
            self.payload_mass += adr.collectedMass+adr.dry_mass+adr.fuel_mass
            
            if self.adr_level_event.triggered:
                self.adr_level_event = env.event()
                self.adrtargetlevel = 0
            
            #print(f"Warehouse {self.id} received {num_satellites} satellite(s) at time {self.env.now}")
            #return num_satellites
            
            matching_sat_ids= [sat for sat in satellites if sat.id in IDs]
            
            # Change the status of collected satellites in the satellites list
            for sat in matching_sat_ids:
                        #print(f"sat {sat.id} has status correctly updated")
                    sat.status += ' warehouse collected'
                    #print(f"sat {sat.id} has status {sat.status}")
                        #sat.status_event.succeed()
                        #sat.status_event = env.event()  # Reset the event for future triggers
                    sat.plane_location = 0
                    sat.warehouse_location = self.id
            
            #change the ADR properties to show that it has offloaded its satellites
            adr.collected_satellites.clear()
            adr.number_collected_satellites = 0
            adr.collectedMass = 0
            
            #now refuel the ADR (if there is enough fuel) and the ADR will become available again
            if self.fuel_mass > adr.fuel_mass_full:
                yield env.timeout(adr.refuel_time)
                adr.status = 'available'
                if not adr.availability_event.triggered:
                    adr.availability_event.succeed()
                adr.fuel_mass = adr.fuel_mass_full
                #print(f"ADR {adr.id} is available again at time {env.now} with fuel mass {adr.fuel_mass_full}")
            else:
                #if there's no fuel for an ADR refill, wait until the warehouse has a resupply event
                self.low_fuel_event.succeed()
                yield self.resupply_event
                yield self.env.timeout(adr.refuel_time)
                self.status = 'available'
                if not adr.availability_event.triggered:
                    self.availability_event.succeed()
                self.fuel_mass = adr.fuel_mass_full
        else:
            #print(f"Warehouse {self.id} does not have enough space to receive {num_satellites} satellites at time {self.env.now}")
            self.adrtargetlevel = target_level
            env.process(self.monitor_collection(env, target_level, adr, IDs,satellites))
            #return 0
        
    def monitor_collection(self, env, target_level, adr, IDs, satellites):
        """if an ADR arrives at a warehouse and there isn't room to deposit its collected satellites, the ADR class function go_home will trigger the warehouse class function monitor_collection.
        When space becomes available, the ADR can offload its collection"""
        while True:
            yield self.adr_level_event
            if self.container.level <= target_level:
                #print(f"monitor collection conditions are satisfied so warehouse {self.id} will not get the satellites with receive_sat")
                self.env.process(self.receive_satellite(env,adr, adr.number_collected_satellites, satellites,IDs,target_level))
            #reset the eventenv, adr,num_satellites, satellites,IDs,target_level
            #self.adr_level_event = env.event()    
            # Change the status of collected satellites in the satellites list
                #for sat_id in IDs:
                #    for sat in satellites:
                #        if sat.id == sat_id:
                #            sat.status += ' warehouse collected'
                            #sat.status_event.succeed()
                            #sat.status_event = env.event()  # Reset the event for future triggers
                        
                #adr.resource = simpy.PriorityResource(env, capacity=adr.max_cap)  # Reset the resource
                #adr.collected_satellites.clear()
                #adr.number_collected_satellites = 0
                #adr.altitude = self.altitude
                #adr.RAAN_0 = self.RAAN_0
                #adr.theta_0 = self.theta_0
        
                #yield env.timeout(refuel_time)
                #adr.status = 'available'
                #adr.fuel_mass = adr.fuel_mass_full
                #print(f"ADR {adr.id} is available again at time {env.now}")
            #yield self.env.timeout(86400)  # Check the container level every 1 day

    def warehouseResupply(self,env,customer,servicer,launcher,satellites,fuel_mass_wh,satellite_drymass,satellite_fuel,electric_fuel_price_per_kg,chemical_fuel_price_per_kg,drift_fuel_per_second,adr_operation_cost,altitude_customer,altitude_parking,sat_dry_mass,sat_ISP,numSatellitesPerPlane,Nparking,ADR_dry_mass,ADR_fuel_mass): 
        """this class process function waits for the warehouse object to reach a critical number of spares, at which point it requests a resupply (if it's not set to de-orbit)
        warehouseResupply also includes the mass and costs necessary to upgrade the warehouse to refuel or repair satellites, if either of those have been requested.
        This function also distinguishes between resupply the warehouse with refurbished satellites and brand new satellites.
        If the customer decides to refurbish previously delivered old satellites, this is deducted from the number of new satellites they must make.
        The customer can decide to make the new satellites refuelable and decide to make the refurbished satellites refuelable."""
        while True:
            #print(f"spare event is {self.below_critical_spare_event}")
            yield self.below_critical_spare_event or self.low_fuel_event
            #yield ((self.container.level-self.collectedSatellites) <= self.criticalSpareLevel)
            if self.deorbit == 0:
                
                #note that the container is full of both ready-to-go spares and old, collected satellites. It requests whatever empty space it has.
                restockNeed = self.max_capacity - self.container.level
                #it also requests a fuel top-off based off of how much it has when full and how much it has currently
                fuelNeed = fuel_mass_wh - self.fuel_mass
                
                
                #if warehouse if refueling satellites, it needs to resupply its satellite fuel
                if self.refuelUpgraded == 1:
                    fuelNeed_satellites = self.refuelUpgradePayload - self.fuel_for_sats
                else:
                    fuelNeed_satellites = 0
                    
                #the total payload includes both new satellites and more fuel
                payload = (restockNeed*(satellite_drymass+satellite_fuel))+fuelNeed+fuelNeed_satellites
                
                #if the launch vehicle has more payload capacity, it is represented with remaining capacity
                remainingCapacity = launcher.maxCapacity - payload
                refuelUp = 0
                repairUp = 0
                
                #after fuel and new satellties have been accounted for, determine if there's space for any wh upgrades, if requested
                if self.upgrade2refuelsats == 1 and remainingCapacity >= self.refuelUpgradePayload:
                    payload += self.refuelUpgradePayload #this needs to include the mass and cost of fuel for refueling satellites (electric propulsion)
                    remainingCapacity -= self.refuelUpgradePayload
                    refuelUp = 1
                if self.upgrade2repairsats == 1 and remainingCapacity >= self.repairUpgradePayload:
                    payload += self.repairUpgradePayload
                    remainingCapacity -= self.repairUpgradePayload
                    repairUp = 1
                    
                #after considering these upgrades, determine if there is room for a new ADR, if one is requested:
                if servicer.newADR == 1 and remainingCapacity >= (ADR_dry_mass+ADR_fuel_mass):
                    satsUnserviced = [sat for sat in satellites if (sat.status == 'retired' or sat.status == 'failed') and sat.assignment == 'unassigned' and sat.deOrbited == 0]
                    #satsUnservicedRAAN = [sat.RAAN_0 for sat in satsUnserviced]

                    # Dictionary to count satellites assigned to each warehouse
                    warehouse_satellite_rev_sum = 0

                    # For each satellite, find the closest warehouse RAAN that does not exceed the satellite RAAN
                    for sat in satsUnserviced:

                        if self.RAAN_0 <= sat.RAAN_0 and ((sat.RAAN_0-self.RAAN_0) < (360/Nparking)):
                            # Increment the count for that warehouse
                            warehouse_satellite_rev_sum += sat.servicePrice
                    
                    """DECISION RULE"""
                    if warehouse_satellite_rev_sum > (sum(sat.servicePrice for sat in satsUnserviced)/Nparking):
                        payload += (ADR_dry_mass+ADR_fuel_mass)
                        remainingCapacity -= (ADR_dry_mass+ADR_fuel_mass)
                        


                
                """When old satellites, collected within warehouses, are returned to Earth, the customer adds them to a list of objects of delivered sats
                if enough time has passed that they were able and willing to refurbish these sats, they are deducted from the number of new satellites they must ship to resupply the warehouse
                These satellites can be made to be refuelable or not (reflected in their refurbishmentPrice)"""
                collectedSatsToGo = 0
                for sat in customer.collectedSatsOnEarth:
                    if (env.now - sat.EarthDeliveryTime) < sat.refurbishmentTime:
                        collectedSatsToGo += 1
                            
                firstSat = satellites[0] #just to extract properties that I need           
                    
                    #DOES CUSTOMER PAY FOR KEEPING THEIR SATS IN THE WAREHOUSE? BY PAYING FOR FUEL PERHAPS?
                    
                """customer takes on cost of new and refurbished satellites. Servicer pays for the warehouse fuel.
                customer pays for the launch costs (compare launch one-by-one vs. batch discount)
                resupply time is taken from the launcher object's property launchTime, which could become an uncertain variable"""
                customer.costs += ((restockNeed-collectedSatsToGo)*firstSat.fullPrice) + collectedSatsToGo*firstSat.refurbishedPrice #just price of the satellites
                servicer.costs += fuelNeed*chemical_fuel_price_per_kg + fuelNeed_satellites*electric_fuel_price_per_kg
                customer.costs += (payload+fuelNeed+fuelNeed_satellites)*launcher.launchCost
                servicer.revenue += (payload*launcher.launchCost)*0.25 #approximates the batch launch benefit -- approximately half as expensive (Jakob), so likely to assume servicer and customer would have profit sharing
                customer.costs += (payload*launcher.launchCost)*0.25
                resupplyTime = launcher.launchTime
                
                    
                """if there was space and a request for refuel and/or repair upgrades, these are reflected in the servicer's costs and added to the resupply time."""
                if self.upgrade2refuelsats == 1 and refuelUp == 1: #if it's requested and there is space on this launch
                    servicer.costs += self.upgrade2refuelCost #(perhaps, upgrade 2 refuel is uncertain based on TRL? )
                    resupplyTime += self.upgrade2refuelTime
                        #self.upgrade2refuelsats = 0
                if self.upgrade2repairsats == 1 and repairUp == 1:
                    servicer.costs += self.upgrade2repairCost #(perhaps, upgrade 2 refuel is uncertain based on TRL? )
                    resupplyTime += self.upgrade2repairTime
                        #self.upgrade2repairsats = 0
                    
                        #trigger launch service
                #print(f"warehouse {self.id} will receive {restockNeed} from launch, which will arrive at time {resupplyTime+self.env.now}")
                self.criticalSpares = False
                self.resupplyOnTheWay = True
                yield self.env.timeout(resupplyTime)
                #print(f"warehouse {self.id} received {restockNeed}. Original container was {self.container.level}")
                
                self.container.put(restockNeed)
                self.numResupply += 1
                if self.upgrade2refuelsats == 1 and refuelUp == 1: #if it's requested and there is space on this launch
                    self.fuel_for_sats = self.refuelUpgradePayload
                #print(f"warehouse {self.id} present container is {self.container.level}")           
                if self.refuelUpgraded == 1:
                    self.fuel_for_sats = self.refuelUpgradePayload
                self.NoSpares = False
                self.criticalSpares = False
                self.resupplyOnTheWay = False
                self.noSpare_event = self.env.event()
                self.below_critical_spare_event = self.env.event()
                self.fuel_mass += fuelNeed #fuel need includes wh fuel
                self.payload_mass += payload #payload includes new satellites, and upgrade mass
                self.resupply_event.succeed
                self.resupply_event = self.env.event()
                #STILL NEEDS LOGIC TO TRACK WHICH OF THE SATELLITES IN THE WAREHOUSE ARE NOW REFUELABLE/REPAIRABLE
                #self.refuelables = 
                #self.repairables =
                
                   
            else: #if the warehouse needs to de-orbit itself
                print(f"warehouse {self.id} is set to be de-orbited")
                
    def warehouse_fullOfOld(self,env):
        while True:
            yield self.fullOfOld_event
            
            
    def warehouseDeorbit(self, env):   
        while True:
            yield self.fullOfOld_event