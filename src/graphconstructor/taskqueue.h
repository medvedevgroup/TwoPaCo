#ifndef _TASK_QUEUE_H_
#define _TASK_QUEUE_H_
#define __STDC_LIMIT_MACROS

#include <mutex>
#include <deque>
#include "common.h"

namespace TwoPaCo
{
	class TaskQueue
	{
	public:
		TaskQueue()
		{

		}

		bool try_push(const Task & task)
		{
			bool ret = false;
			mutex_.lock();

			if (queue_.size() < capacity_)
			{
				queue_.push_back(task);
				ret = true;
			}

			mutex_.unlock();
			return ret;
		}

		bool try_pop(Task& task)
		{
			bool ret = false;
			mutex_.lock();
			if (queue_.size() > 0)
			{
				task = queue_.front();
				queue_.pop_front();
				ret = true;
			}

			mutex_.unlock();
			return ret;
		}

		void set_capacity(size_t capacity)
		{
			mutex_.lock();
			capacity_ = capacity;
			mutex_.unlock();
		}

		size_t size()
		{
			mutex_.lock();
			size_t ret = queue_.size();
			mutex_.unlock();
			return ret;
		}

		size_t capacity() const
		{
			return capacity_;
		}

	private:
		size_t capacity_;
		std::mutex mutex_;
		std::deque<Task> queue_;
		DISALLOW_COPY_AND_ASSIGN(TaskQueue);
	};

	typedef std::unique_ptr<TaskQueue> TaskQueuePtr;
}

#endif
