// ch13_arduino_laser.ino
//
// 原书 laser.ino：经 TB6612FNG 双电机驱动板驱动两台直流电机（各带一面旋转镜反射激光）
// 引脚映射：STBY=10；电机#1（A）PWMA=3、AIN1=9、AIN2=8；电机#2（B）PWMB=5、BIN1=11、BIN2=12
// 串口协议：'H' 头校验 + 4 字节（speed1, dir1, speed2, dir2），9600 波特

int STBY = 10; //standby

// Motor #1
int PWMA = 3; //Speed control
int AIN1 = 9; //Direction
int AIN2 = 8; //Direction

// Motor #2
int PWMB = 5; //Speed control
int BIN1 = 11; //Direction
int BIN2 = 12; //Direction

void setup(){

  pinMode(STBY, OUTPUT);

  pinMode(PWMA, OUTPUT);
  pinMode(AIN1, OUTPUT);
  pinMode(AIN2, OUTPUT);

  pinMode(PWMB, OUTPUT);
  pinMode(BIN1, OUTPUT);
  pinMode(BIN2, OUTPUT);

  // initialize the serial communication:
  Serial.begin(9600);
}

// main loop that reads motor data sent by laser.py
void loop()
{
  // The data sent is of the form:
  // 'H' (header), speed1, dir1, speed2, dir2
  if (Serial.available() >= 5) {
    if (Serial.read() == 'H') {
      // read next 4 bytes
      byte s1 = Serial.read();
      byte d1 = Serial.read();
      byte s2 = Serial.read();
      byte d2 = Serial.read();

      // stop motor if both speeds are 0
      if (s1 == 0 && s2 == 0) {
        stop();
      }
      else {
        // set speed and direction of both motors
        move(0, s1, d1);
        move(1, s2, d2);
      }
      // small pause
      delay(20);
    }
    else {
      // invalid data, stop motors
      stop();
    }
  }
  else {
    // no data
    delay(250);
  }
}

// set motor speed and direction
// motor: #1 -> 1, #2 -> 0
// direction: 1/0
void move(int motor, int speed, int direction)
{
  // disable standby
  digitalWrite(STBY, HIGH);

  boolean inPin1 = LOW;
  boolean inPin2 = HIGH;

  if (direction == 1) {
    inPin1 = HIGH;
    inPin2 = LOW;
  }

  if (motor == 1) {
    digitalWrite(AIN1, inPin1);
    digitalWrite(AIN2, inPin2);
    analogWrite(PWMA, speed);
  }
  else {
    digitalWrite(BIN1, inPin1);
    digitalWrite(BIN2, inPin2);
    analogWrite(PWMB, speed);
  }
}

// bring motors to a stop
void stop(){
  // enable standby
  digitalWrite(STBY, LOW);
}
