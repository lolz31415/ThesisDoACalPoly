// https://forum.arduino.cc/t/433mhz-transmitter-on-esp32/1115389/4

// ESP32 433MHz receiver test 1

#include <RH_ASK.h>  // Include RadioHead Amplitude Shift Keying Library
#include <SPI.h>     // Include dependant SPI Library

// speed (bits per second), rx_pin, tx_pin
// GPIO pin_22 -> D22 (2nd to top right side of ESP-Wroom23)
RH_ASK rf_driver(2000, 21, 23);   // ESP32 Create Amplitude Shift Keying Object
// can set to promiscuous mode:
// Tells the receiver to accept messages with any TO address, not just messages addressed to thisAddress or the broadcast address

void setup() {
  Serial.begin(115200);
  //rf_driver.init();
  
  delay(4000);
  Serial.println("ESP32 433MHz receiver");
  if (RH_PLATFORM == RH_PLATFORM_ESP32)
    Serial.println("RH_PLATFORM_ESP32");
  delay(5000);
  Serial.println("Transmitter: rf_driver initialising");
  if (!rf_driver.init()) {
    Serial.println("init failed");
    while (1) delay(1000);
  }
  Serial.println("Transmitter: rf_driver initialised");
  
}

// transmit packet every 5 seconds
void loop() {
  Serial.println("Transmitting packet");
  const char *msg = "Hello World";
  rf_driver.send((uint8_t *)msg, strlen(msg)+1);
  rf_driver.waitPacketSent();
  delay(1000);
}